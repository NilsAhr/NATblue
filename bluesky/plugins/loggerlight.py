""" BlueSky loggerlight - lightweight flight status logger
Reduced-column version of loggerff for smaller log files.
Created on 2025
@author: ahre_ni, Nils Ahrenhold
"""
import numpy as np
from bluesky import core, stack, traf, sim, navdb
from bluesky.core import Entity
from bluesky.tools import datalog, geo
from bluesky.tools.aero import ft, kts, nm, fpm, vtas2cas
from bluesky.tools.position import txt2pos
import bluesky as bs

flstheader = \
    'simt,' + \
    'flightid,' + \
    'ac_type,' + \
    'flighttime,' + \
    'distanceflown,' + \
    'totalfuel,' + \
    'latitude,' + \
    'longitude,' + \
    'positivefuelflow,' + \
    'rawfuelflow,' + \
    'thrust_effective,' + \
    'altitude,' + \
    'tas,' + \
    'gs,' + \
    'cas,' + \
    'mach,' + \
    'vs,' + \
    'n_active_conflicts,' + \
    'n_active_intrusions,' + \
    'flight_phase' + '\n'


### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
def init_plugin():
    ''' Plugin initialisation function. '''
    global loggerlight
    loggerlight = Loggerlight()

    config = {
        'plugin_name':     'LOGGERLIGHT',
        'plugin_type':     'sim',
    }

    stackfunctions = {
        'STARTLOGLIGHT': [
            'STARTLOGLIGHT',
            '',
            loggerlight.start_log,
            'Starts the lightweight flight status logger'
        ]
    }

    return config, stackfunctions


class Loggerlight(Entity):
    ''' Lightweight flight status logger. '''
    def __init__(self):
        super().__init__()

        self.d = 10  # distance parameter [nm] for deleting aircraft automatically
        self.sim_name = stack.get_scenname()

        # The FLST LOG
        self.flst = datalog.crelog('FLSTLOG_LOGGERLIGHT', None, flstheader)

        with self.settrafarrays():
            self.distance2D = np.array([])
            self.create_time = np.array([])
            self.total_fuel = np.array([])
            self.last_update_time = np.array([])

    def reset(self):
        super().reset()
        self.sim_name = None

    def create(self, n=1):
        ''' Called automatically when new aircraft are created. '''
        super().create(n)
        self.create_time[-n:] = sim.simt
        self.total_fuel[-n:] = 0.0
        self.last_update_time[-n:] = sim.simt

    @core.timed_function(name='LOGGERLIGHT', dt=1.0)
    def update(self, dt):
        """
        1. Handle aircraft deletion
        2. Early exit if no aircraft
        3. Update tracking arrays
        4. Calculate fuel
        5. Log flight status
        """
        n_current = len(traf.id)

        #########################################################
        ################## 1. AIRCRAFT DELETION ################
        #########################################################
        if traf.ntraf > 0:
            lat_ac = np.array(traf.lat)
            lon_ac = np.array(traf.lon)
            dest_names = np.array(traf.ap.dest)

            lat_dest = np.zeros(traf.ntraf)
            lon_dest = np.zeros(traf.ntraf)
            valid_dest = np.ones(traf.ntraf, dtype=bool)

            for idx, apname in enumerate(dest_names):
                apidx = bs.navdb.getaptidx(apname)

                if apidx < 0:
                    if bs.traf.ap.route[idx].nwp > 0:
                        lat_dest[idx] = bs.traf.ap.route[idx].wplat[-1]
                        lon_dest[idx] = bs.traf.ap.route[idx].wplon[-1]
                    else:
                        lat_dest[idx] = traf.lat[idx]
                        lon_dest[idx] = traf.lon[idx]

                    success, posobj = txt2pos(apname, lat_dest[idx], lon_dest[idx])
                    if success:
                        lat_dest[idx] = posobj.lat
                        lon_dest[idx] = posobj.lon
                    else:
                        valid_dest[idx] = False
                else:
                    lat_dest[idx] = bs.navdb.aptlat[apidx]
                    lon_dest[idx] = bs.navdb.aptlon[apidx]

            _, distances = geo.qdrdist(
                lat_ac[valid_dest], lon_ac[valid_dest],
                lat_dest[valid_dest], lon_dest[valid_dest])

            to_delete_valid = np.where(distances < self.d)[0]
            to_delete = np.nonzero(valid_dest)[0][to_delete_valid]

            for idx in to_delete:
                cs = traf.id[idx]
                traf.delete(idx)
                ac_act = traf.ntraf
                print(f"LOGGERLIGHT - {self.sim_name}: {cs} landed at {sim.simt}; active aircraft: {ac_act}")

        #########################################################
        ################## 2. EARLY EXIT CHECK ################
        #########################################################
        if sim.simt >= 24 * 3600:
            if traf.ntraf == 0:
                name = self.sim_name or "NAME NOT DEFINED"
                print(f"END of simulation: {name} at: {sim.simt}s = {sim.simt/3600}h")
                stack.stack('reset')
                return

        if traf.ntraf == 0:
            return

        #########################################################
        ################## 3. UPDATE TRACKING ARRAYS ##########
        #########################################################
        n_current = traf.ntraf

        if n_current > 0:
            self.distance2D[:n_current] += dt * traf.gs[:n_current]

        #########################################################
        ################## 4. FUEL CALCULATIONS ###############
        #########################################################
        try:
            raw_fuelflow = traf.perf.fuelflow[:n_current]
            positive_fuelflow = np.maximum(0, raw_fuelflow)

            # Use thrust_effective (matches fuel flow) for logging
            thrust_eff = np.zeros(n_current)
            if hasattr(traf.perf, 'thrust_effective') and len(traf.perf.thrust_effective) >= n_current:
                thrust_eff = traf.perf.thrust_effective[:n_current]
            elif hasattr(traf.perf, 'thrust') and len(traf.perf.thrust) >= n_current:
                # Fallback to raw thrust if thrust_effective not available
                thrust_eff = traf.perf.thrust[:n_current]

            fuel_consumed_this_step = positive_fuelflow * 1.0
            self.total_fuel[:n_current] += fuel_consumed_this_step
            self.last_update_time[:n_current] = sim.simt

        except (AttributeError, IndexError):
            raw_fuelflow = np.full(n_current, np.nan)
            positive_fuelflow = np.full(n_current, 0.0)
            thrust_eff = np.full(n_current, np.nan)

        #########################################################
        ################## 5. FLIGHT STATUS LOGGING ###########
        #########################################################
        n_active_conflicts = len(traf.cd.confpairs_unique)
        n_active_intrusions = len(traf.cd.lospairs_unique)

        flight_phase = np.zeros(n_current)
        if hasattr(traf.perf, 'phase') and len(traf.perf.phase) >= n_current:
            flight_phase = traf.perf.phase[:n_current]

        if n_current > 0:
            self.flst.log(
                traf.id[:n_current],                               # flightid [-]
                traf.type[:n_current],                             # ac_type [-]
                (sim.simt - self.create_time[:n_current]),         # flighttime [s]
                (traf.distflown[:n_current] / nm),                 # distanceflown [nm]
                self.total_fuel[:n_current],                       # totalfuel [kg]
                traf.lat[:n_current],                              # latitude [deg]
                traf.lon[:n_current],                              # longitude [deg]
                positive_fuelflow,                                 # positivefuelflow [kg/s]
                raw_fuelflow,                                      # rawfuelflow [kg/s]
                thrust_eff,                                        # thrust_effective [N]
                (traf.alt[:n_current] / ft),                       # altitude [ft]
                (traf.tas[:n_current] / kts),                      # tas [kts]
                (traf.gs[:n_current] / kts),                       # gs [kts]
                (traf.cas[:n_current] / kts),                      # cas [kts]
                traf.M[:n_current],                                # mach [-]
                (traf.vs[:n_current] / fpm),                       # vs [fpm]
                n_active_conflicts,                                # n_active_conflicts [-]
                n_active_intrusions,                               # n_active_intrusions [-]
                flight_phase.astype(int),                          # flight_phase [1-6]
            )

    def start_log(self):
        print(f"LOGGERLIGHT: flight status logger started: {self.sim_name}")
        self.flst.start()
        return True, f'LOGGERLIGHT is ON for simulation: {self.sim_name}.'
