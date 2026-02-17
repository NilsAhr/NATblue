""" BlueSky logger_ff for conflict and flight state in North Atlantic
Created on 2025 Sep
@author: ahre_ni, Nils Ahrenhold
"""
from random import randint
import numpy as np
# Import the global bluesky objects. Uncomment the ones you need
from bluesky import core, stack, traf, sim, navdb  #settings, scr, tools
from bluesky.core import Entity
from bluesky.tools import datalog, geo
from bluesky.tools.aero import ft,kts,nm,fpm, vtas2cas, vtas2mach
from bluesky.tools.position import txt2pos
import bluesky as bs
from bluesky.traffic.performance.perfbase import PerfBase

flstheader = \
    'simt,' + \
    'flightid,' + \
    'ac_type,' + \
    'flighttime,' + \
    'currentmass,' + \
    'distanceflown,' + \
    'totalfuel,' + \
    'latitude,' + \
    'longitude,' + \
    'spawntime,' + \
    'actualdistance2D,' + \
    'actualdistance3D,' + \
    'workdone,' + \
    'positivefuelflow,' + \
    'rawfuelflow,' + \
    'thrust,' + \
    'altitude,' + \
    'tas,' + \
    'gs,' + \
    'cas,' + \
    'mach,' + \
    'vs,' + \
    'heading,' + \
    'asasactive,' + \
    'pilotalt,' + \
    'pilottas,' + \
    'pilothdg,' + \
    'pilotvs,' + \
    'n_active_conflicts,' + \
    'n_active_intrusions,' + \
    'gsnorth,' + \
    'gseast,' + \
    'windnorth,' + \
    'windeast,' + \
    'headwind,' + \
    'crosswind,' + \
    'pilotcas,' + \
    'pilotmach,' + \
    'selalt,' + \
    'selvs,' + \
    'swlnav,' + \
    'swvnav,' + \
    'swvnavspd,' + \
    'swats,' + \
    'throttle,' + \
    'temp,' + \
    'rho,' + \
    'flight_phase,' + \
    'drag' + '\n'

confheader = \
    'simt[s],' + \
    'ac1,' + \
    'ac2,' + \
    'latitude_ac1[deg],' + \
    'longitude_ac1[deg],' + \
    'altitude_ac1[ft],' + \
    'latitude_ac2[deg],' + \
    'longitude_ac2[deg],' + \
    'altitude_ac2[ft],' + \
    'heading_ac1[deg],' + \
    'heading_ac2[deg],' + \
    'vs_ac1[fpm],' + \
    'vs_ac2[fpm],' + \
    'dcpa[nm],' + \
    'tcpa[sec],' + \
    'tLOS[sec],' + \
    'qdr[deg],' + \
    'dist[nm],' +\
    'dalt_min[ft]' + \
    'tinconf[sec],' + \
    'toutconf[sec],' + \
    'duration[s],' + \
    'intrusion' + '\n'


### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
def init_plugin():
    ''' Plugin initialisation function. '''
    # Instantiate our example entity
    global loggerff
    loggerff = Loggerff()

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'LOGGERFF',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim',
        }
    
    stackfunctions = {
        'STARTLOG': [
            'STARTLOG',
            '',
            loggerff.start_log,
            'Starts the flight status and conflict loggerff'
        ],
        'SETMASS': [
            'SETMASS acid,mass',
            '[acid,float]',
            loggerff.setmass,
            'Set aircraft mass [kg] after creation, e.g. SETMASS callsign, 257656.4'
        ]
        }

    # init_plugin() should always return a configuration dict.
    return config, stackfunctions


class Loggerff(Entity):
    ''' Example new entity object for BlueSky. '''
    def __init__(self):
        super().__init__()

        # Parameters for conflict count and statistics
        self.duration = {}              #dict for duration values
        self.d = 10                     #distance parameter for deleting aircraft automatically 

        self.mass_warning_flightids = set()  # Track which aircraft IDs we've warned about

        # basic conflict parameters
        self.init_lat1 = {}
        self.init_lon1 = {}
        self.init_alt1 = {}
        self.init_hdg1 = {}
        self.init_vs1 = {}
        self.init_lat2 = {}
        self.init_lon2 = {}
        self.init_alt2 = {}
        self.init_hdg2 = {}
        self.init_vs2 = {}

        # check for intrusion
        self.intrusion_occurred = {}    # Track if an intrusion occurred for each conflict

        # new severity parameters
        self.dcpa = {}    # Track minimum dcpa for each conflict
        self.dalt = {}    # Track minimum dalt for each conflict
        self.tLOS = {}
        self.dist = {}
        self.qdr = {}
        self.tcpa = {}
        self.tinconf = {}     # Track conflict start times for each conflict
        self.toutconf = {}    # Track conflict end times for each conflict
        self.sim_name = stack.get_scenname()

        # Performance model access
        #self.perf     = PerfBase()

        # The FLST & CONF LOGGERFF
        self.flst = datalog.crelog('FLSTLOG_LOGGERFF', None, flstheader)
        self.conflog = datalog.crelog('CONFLOG_LOGGERFF', None, confheader)

        with self.settrafarrays():
            self.distance2D = np.array([])
            self.distance3D = np.array([])
            self.create_time = np.array([])
            self.total_fuel = np.array([]) # total fuel [kg]
            self.last_update_time = np.array([]) # last update time [s] for integration

    def reset(self):
        super().reset()
        self.duration = {}

        self.mass_warning_flightids = set()  # Reset mass warning tracking

        #basic conflict parameters
        self.init_lat1 = {}
        self.init_lon1 = {}
        self.init_alt1 = {}
        self.init_hdg1 = {}
        self.init_vs1 = {}
        self.init_lat2 = {}
        self.init_lon2 = {}
        self.init_alt2 = {}
        self.init_hdg2 = {}
        self.init_vs2 = {}

        self.intrusion_occurred = {}

        # reset performance model
        #self.perf.reset()

        # severity parameters
        self.dcpa = {}
        self.dalt = {}
        self.tLOS = {}
        self.dist = {}
        self.qdr = {}
        self.tcpa = {}
        self.tinconf = {}
        self.toutconf = {}
        self.sim_name = None


    def create(self, n=1):
        ''' This function gets called automatically when new aircraft are created. '''
        super().create(n)        
        self.create_time[-n:] = sim.simt
        # Initialize fuel tracking for new aircraft
        self.total_fuel[-n:] = 0.0
        self.last_update_time[-n:] = sim.simt

    def delete(self, idx):
        # Clean up when aircraft are deleted (called automatically by BlueSky)
        super().delete(idx)
        
        # Clean up mass warning tracking for deleted aircraft
        if idx < len(traf.id):
            deleted_flightid = traf.id[idx]
            self.mass_warning_flightids.discard(deleted_flightid)
            
            # Optional: Log cleanup for debugging
            # print(f"LOGGERFF: Cleaned up tracking for deleted aircraft: {deleted_flightid}")

    @core.timed_function(name='LOGGERFF', dt=1.0)
    def update(self, dt):
        """
        1. Handle aircraft deletion first
        2. Early exit if no aircraft
        3. Update tracking arrays for remaining aircraft  
        4. Calculate derived values
        5. Log everything with consistent array sizes
        """
        n_current = len(traf.id)
                
        #########################################################
        ################## 1. AIRCRAFT DELETION ################
        #########################################################
        # Do this FIRST before any calculations
        
        if traf.ntraf > 0:  # Only if there are aircraft
            # Extract the latitudes and longitudes of all aircraft
            lat_ac = np.array(traf.lat)       # All aircraft latitudes
            lon_ac = np.array(traf.lon)       # All aircraft longitudes
            dest_names = np.array(traf.ap.dest)  # Destination airport names for each aircraft

            # Initialize destination arrays
            lat_dest = np.zeros(traf.ntraf)
            lon_dest = np.zeros(traf.ntraf)

            # Boolean array to track which aircraft have valid destinations
            valid_dest = np.ones(traf.ntraf, dtype=bool)  # Start by assuming all destinations are valid
            
            # Precompute destination latitudes and longitudes for all aircraft
            for idx, apname in enumerate(dest_names):
                apidx = bs.navdb.getaptidx(apname)
                
                if apidx < 0:
                    # No valid destination, use the last waypoint or current position
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
                        valid_dest[idx] = False  # Mark destination as invalid
                else:
                    # Valid destination airport
                    lat_dest[idx] = bs.navdb.aptlat[apidx]
                    lon_dest[idx] = bs.navdb.aptlon[apidx]

            # Now compute the distances between each aircraft and its destination
            _, distances = geo.qdrdist(lat_ac[valid_dest], lon_ac[valid_dest], lat_dest[valid_dest], lon_dest[valid_dest])

            # Vectorized deletion: Find indices where the distance is less than 10nm
            to_delete_valid = np.where(distances < self.d)[0]

            # Map valid indices back to the original indices (since we filtered by valid_dest)
            to_delete = np.nonzero(valid_dest)[0][to_delete_valid]

            # Delete the aircraft all at once
            for idx in to_delete:
                cs = traf.id[idx]
                traf.delete(idx)
                ac_act = traf.ntraf
                print(f"FLST LOGGERFF - {self.sim_name}: {cs} landed at {sim.simt}; active aircraft: {ac_act}")

        #########################################################
        ################## 2. EARLY EXIT CHECK ################
        #########################################################
        
        # Check if simulation should end
        if sim.simt >= 48 * 3600:  # 24 hours in seconds
            if traf.ntraf == 0:
                if self.sim_name is None:
                    print(f"END of simulation: NAME NOT DEFINED at: {sim.simt}seconds = {sim.simt/3600}hours")
                else:
                    print(f"END of simulation: {self.sim_name} at: {sim.simt}seconds = {sim.simt/3600}hours")
                stack.stack('reset')
                return  # EXIT EARLY - no aircraft to process
        
        # If no aircraft, don't do any calculations
        if traf.ntraf == 0:
            return

        #########################################################
        ################## 3. UPDATE TRACKING ARRAYS ##########
        #########################################################
        
        # Ensure all our tracking arrays match current aircraft count
        n_current = traf.ntraf
        
        # Calculate time since last update for each aircraft
        dt_array = 1.0

        # Update distance tracking
        if n_current > 0:
            resultantspd = np.sqrt(traf.gs[:n_current] ** 2 + traf.vs[:n_current] ** 2)
            self.distance2D[:n_current] += dt * traf.gs[:n_current]
            self.distance3D[:n_current] += dt * resultantspd

        #########################################################
        ################## 4. FUEL CALCULATIONS ###############
        #########################################################
        
        # Get fuel flow data (ensure correct size)
        try:
            raw_fuelflow = traf.perf.fuelflow[:n_current]
            positive_fuelflow = np.maximum(0, raw_fuelflow)

            # Thrust - use thrust_effective (matches fuel flow) if available
            thrust = np.zeros(n_current)
            if hasattr(traf.perf, 'thrust_effective') and len(traf.perf.thrust_effective) >= n_current:
                thrust = traf.perf.thrust_effective[:n_current]
            elif hasattr(traf.perf, 'thrust') and len(traf.perf.thrust) >= n_current:
                thrust = traf.perf.thrust[:n_current]
            
            # Add this debug code for mass
            current_mass = np.zeros(n_current)
            if hasattr(traf.perf, 'mass') and len(traf.perf.mass) >= n_current:
                current_mass = traf.perf.mass[:n_current]
                #print(f"Time: {sim.simt:.1f}, Mass in logger: {current_mass[0]:.2f}, Fuelflow: {raw_fuelflow[0]:.5f}")

            # Integrate fuel consumption
            fuel_consumed_this_step = positive_fuelflow * dt_array
            self.total_fuel[:n_current] += fuel_consumed_this_step
            
            # Update timestamp
            self.last_update_time[:n_current] = sim.simt
            
        except (AttributeError, IndexError):
            # Fallback if fuel flow not available
            raw_fuelflow = np.full(n_current, np.nan)
            positive_fuelflow = np.full(n_current, 0.0)
            thrust = np.full(n_current, np.nan)
        
        
        #########################################################
        ################## 5. CONFLICT PROCESSING #############
        #########################################################        
        
        # Always update min values for all ongoing conflicts
        idxdict = {frozenset(v): i for i, v in enumerate(traf.cd.confpairs)}
        for confpair_frozen in traf.cd.confpairs_unique:
            if confpair_frozen in idxdict:
                i = idxdict[confpair_frozen]
                dcpa_now = np.asarray(traf.cd.dcpa)[i]
                dalt_now = np.asarray(traf.cd.dalt)[i]
                tcpa_now = np.asarray(traf.cd.tcpa)[i]
                tLOS_now = np.asarray(traf.cd.tLOS)[i]
                qdr_now = np.asarray(traf.cd.qdr)[i]
                dist_now = np.asarray(traf.cd.dist)[i]
                # Store latest values for each conflict pair
                self.dcpa[confpair_frozen] = dcpa_now
                self.dalt[confpair_frozen] = dalt_now # already updated for minimum in statebased.py
                self.tLOS[confpair_frozen] = tLOS_now
                self.qdr[confpair_frozen] = qdr_now
                self.tcpa[confpair_frozen] = tcpa_now

                # Track minimum distance of conflict pairs
                if confpair_frozen not in self.dist:
                    self.dist[confpair_frozen] = dist_now
                else:
                    self.dist[confpair_frozen] = min(self.dist[confpair_frozen], dist_now)
                
                # Set tinconf only once, when conflict starts
                if confpair_frozen not in self.tinconf:
                    self.tinconf[confpair_frozen] = sim.simt

                # Check if this pair is currently an intrusion
                if confpair_frozen in traf.cd.lospairs_unique:
                    self.intrusion_occurred[confpair_frozen] = True
                elif confpair_frozen not in self.intrusion_occurred:
                    self.intrusion_occurred[confpair_frozen] = False

        ##########################################################
        ############### CONFLICT DURATION & LOGGING ##############
        ########################################################## 
        # 
        # - conflicts still ongoing when simulation ends are not logged
        #       
        # Loop through each conflict pair in traf.cd.confpairs
        for pair in traf.cd.confpairs:
        ########## NEW CONFLICT PAIR DETECTION ##########
        # If the conflict pair is new, initialize its duration to 1
            # Use frozenset to store unique pairs in duration dictionary
            unique_pair = frozenset(pair)
            if unique_pair not in self.duration:
                self.duration[unique_pair] = 1
                # Store initial parameters
                ac1, ac2 = tuple(pair)
                idx1 = traf.id2idx(ac1)
                idx2 = traf.id2idx(ac2)
                self.init_lat1[unique_pair] = traf.lat[idx1]
                self.init_lon1[unique_pair] = traf.lon[idx1]
                self.init_alt1[unique_pair] = traf.alt[idx1]
                self.init_hdg1[unique_pair] = traf.hdg[idx1]
                self.init_vs1[unique_pair] = traf.vs[idx1]/fpm

                self.init_lat2[unique_pair] = traf.lat[idx2]
                self.init_lon2[unique_pair] = traf.lon[idx2]
                self.init_alt2[unique_pair] = traf.alt[idx2]
                self.init_hdg2[unique_pair] = traf.hdg[idx2]
                self.init_vs2[unique_pair] = traf.vs[idx2]/fpm
            else:
                # If the conflict pair already exists, its duration has already started tracking
                continue
        ########## CONFLICT DURATION UPDATE ##########
        # Now, loop through each pair in the duration dictionary
        for pair, duration in list(self.duration.items()):  # Use list() to allow modifying the dict during iteration
            confpair_frozen = frozenset(pair)

            # Check if the pair is still part of traf.cd.confpairs_unique
            if confpair_frozen in traf.cd.confpairs_unique:
                # Conflict is still ongoing, increment the duration by 1
                self.duration[pair] += 1
        ############## CONFLICT END DETECTION ##########
            else:
                # Conflict has ended, log the pair and its duration
                ac1, ac2 = tuple(pair)
                idx1 = traf.id2idx(ac1)
                idx2 = traf.id2idx(ac2)
                lat_1 = self.init_lat1[confpair_frozen]
                lon_1 = self.init_lon1[confpair_frozen]
                alt_1 = self.init_alt1[confpair_frozen]
                hdg_1 = self.init_hdg1[confpair_frozen]
                vs_1  = self.init_vs1[confpair_frozen]

                lat_2 = self.init_lat2[confpair_frozen]
                lon_2 = self.init_lon2[confpair_frozen]
                alt_2 = self.init_alt2[confpair_frozen]
                hdg_2 = self.init_hdg2[confpair_frozen]
                vs_2  = self.init_vs2[confpair_frozen]

        ######## CONF LOGGING #########################
                # Log the parameters only once, at conflict end
                self.conflog.log(
                    ac1, ac2, lat_1, lon_1, alt_1/ft, lat_2, lon_2, alt_2/ft, hdg_1, hdg_2, vs_1, vs_2, #initial parameters
                    self.dcpa[confpair_frozen] / nm, self.tcpa[confpair_frozen], self.tLOS[confpair_frozen], #conflict parameters
                    self.qdr[confpair_frozen], self.dist[confpair_frozen] / nm, self.dalt[confpair_frozen] / ft,
                    self.tinconf[confpair_frozen], sim.simt, self.duration[pair],
                    int(self.intrusion_occurred[confpair_frozen]) #1 if intrusion occurred, 0 otherwise
                )                
                del self.duration[pair]
                self.toutconf[confpair_frozen] = sim.simt
                # Optionally, clean up min_dcpa, min_dalt, tinconf, toutconf for this pair

        
        #########################################################
        ################## 6. FLIGHT STATUS LOGGING ###########
        #########################################################
        
        # Count active conflicts and intrusions (these are scalars, not arrays)
        n_active_conflicts = len(traf.cd.confpairs_unique)
        n_active_intrusions = len(traf.cd.lospairs_unique)

        # Precompute derived debug signals
        # Track and wind to headwind/xwind in m/s
        trk_rad = np.radians(traf.trk[:n_current])
        wn, we = traf.windnorth[:n_current], traf.windeast[:n_current]
        headwind = -(wn * np.cos(trk_rad) + we * np.sin(trk_rad))
        crosswind = (-wn * np.sin(trk_rad) + we * np.cos(trk_rad))

        # Commanded speed regime from pilottas
        pilot_tas = traf.aporasas.tas[:n_current]  # [m/s]
        pilot_cas_kn = vtas2cas(pilot_tas, traf.alt[:n_current]) / kts
        pilot_mach = vtas2mach(pilot_tas, traf.alt[:n_current])

        # Autopilot/ASAS selected targets
        selalt_ft = traf.aporasas.alt[:n_current] / ft
        selvs_fpm = traf.aporasas.vs[:n_current] / fpm

        # Mode flags & throttle
        swlnav = traf.swlnav[:n_current]
        swvnav = traf.swvnav[:n_current]
        swvnavspd = traf.swvnavspd[:n_current]
        swats = traf.swats[:n_current] if len(traf.swats) >= n_current else np.zeros(n_current, dtype=bool)
        throttle = traf.thr[:n_current] if len(traf.thr) >= n_current else np.full(n_current, -999.0)

        # Atmosphere (already maintained in traffic)
        Temp = traf.Temp[:n_current]
        rho = traf.rho[:n_current]

        # Flight phase and drag from performance model (for diagnostics)
        flight_phase = np.zeros(n_current)
        drag = np.zeros(n_current)
        if hasattr(traf.perf, 'phase') and len(traf.perf.phase) >= n_current:
            flight_phase = traf.perf.phase[:n_current]
        if hasattr(traf.perf, 'D') and len(traf.perf.D) >= n_current:
            drag = traf.perf.D[:n_current]

        # Ensure ALL arrays have exactly n_current elements and log
        if n_current > 0:
            self.flst.log(
                traf.id[:n_current],                               # flightid [-]
                traf.type[:n_current],                             # ac_type [-]
                (sim.simt - self.create_time[:n_current]),         # flighttime [s]
                current_mass,                                      # currentmass [kg]
                (traf.distflown[:n_current] / nm),                 # distanceflown [nm]
                self.total_fuel[:n_current],                       # totalfuel [kg]
                traf.lat[:n_current],                              # latitude [deg]
                traf.lon[:n_current],                              # longitude [deg]
                self.create_time[:n_current],                      # spawntime [s]
                self.distance2D[:n_current],                       # actualdistance2D [m]
                self.distance3D[:n_current],                       # actualdistance3D [m]
                (traf.work[:n_current] * 1e-6),                    # workdone [MJ]
                positive_fuelflow,                                 # positivefuelflow [kg/s]
                raw_fuelflow,                                      # rawfuelflow [kg/s]
                thrust,                                            # thrust [N]
                (traf.alt[:n_current] / ft),                       # altitude [ft]
                (traf.tas[:n_current] / kts),                      # tas [kts]
                (traf.gs[:n_current] / kts),                       # gs [kts]
                (traf.cas[:n_current] / kts),                      # cas [kts]
                traf.M[:n_current],                                # mach [-]
                (traf.vs[:n_current] / fpm),                       # vs [fpm]
                traf.hdg[:n_current],                              # heading [deg]
                traf.cr.active[:n_current],                        # asasactive [bool/int]
                (traf.aporasas.alt[:n_current] / ft),              # pilotalt [ft]
                (traf.aporasas.tas[:n_current] / kts),             # pilottas [kts]
                traf.aporasas.hdg[:n_current],                     # pilothdg [deg]
                (traf.aporasas.vs[:n_current] / fpm),              # pilotvs [fpm]
                n_active_conflicts,                                # n_active_conflicts [-]
                n_active_intrusions,                               # n_active_intrusions [-]
                (traf.gsnorth[:n_current] / kts),                  # gsnorth [kts]
                (traf.gseast[:n_current] / kts),                   # gseast [kts]
                (traf.windnorth[:n_current] / kts),                # windnorth [kts]
                (traf.windeast[:n_current] / kts),                 # windeast [kts]
                (headwind / kts),                                  # headwind [kts]
                (crosswind / kts),                                 # crosswind [kts]
                pilot_cas_kn,                                      # pilotcas [kts]
                pilot_mach,                                        # pilotmach [-]
                selalt_ft,                                         # selalt [ft]
                selvs_fpm,                                         # selvs [fpm]
                swlnav.astype(int),                                # swlnav [0/1]
                swvnav.astype(int),                                # swvnav [0/1]
                swvnavspd.astype(int),                             # swvnavspd [0/1]
                swats.astype(int),                                 # swats [0/1]
                throttle,                                          # throttle [-]
                Temp,                                              # temp [K]
                rho,                                               # rho [kg/m³]
                flight_phase.astype(int),                          # flight_phase [1-6]
                drag                                               # drag [N]
            )    
            
    def start_log(self):
        print (f"LOGGERFF: FLST and CONF loggerff started: {self.sim_name}")
        self.flst.start()
        self.conflog.start()
        return True, f'FLST and CONF loggerff is ON for simulation: {self.sim_name}.'

    def setmass(self, idx, mass_kg):
        """Set aircraft mass [kg] in the BADA performance model.
        
        Usage in scenario file:
            SETMASS callsign, mass_in_kg
        
        This sets the initial mass for the aircraft after CRE, overriding
        the default 95% MTOW used by BADA.
        """
        # Validate mass
        if mass_kg <= 0:
            return False, f'SETMASS: Mass must be positive, got {mass_kg} kg.'

        # Get the flight ID for messaging
        acid = traf.id[idx]

        # Check if BADA performance model is active and has mass array
        if not hasattr(traf.perf, 'mass') or len(traf.perf.mass) == 0:
            return False, f'SETMASS {acid}: No performance model with mass data active. Is PERF BADA loaded?'

        # Validate against BADA mass limits if available
        if hasattr(traf.perf, 'mmin') and hasattr(traf.perf, 'mmax'):
            mmin = traf.perf.mmin[idx]
            mmax = traf.perf.mmax[idx]
            if mass_kg < mmin:
                return False, f'SETMASS {acid}: Mass {mass_kg:.1f} kg is below OEW limit ({mmin:.1f} kg).'
            if mass_kg > mmax:
                return False, f'SETMASS {acid}: Mass {mass_kg:.1f} kg exceeds MTOW limit ({mmax:.1f} kg).'

        # Set the mass
        old_mass = traf.perf.mass[idx]
        traf.perf.mass[idx] = mass_kg

        print(f'SETMASS {acid}: Mass set from {old_mass:.1f} kg to {mass_kg:.1f} kg')
        return True, f'SETMASS {acid}: Mass set to {mass_kg:.1f} kg (was {old_mass:.1f} kg).'
