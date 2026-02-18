""" BlueSky loggercompact — minimal flight-state & conflict logger
    for North Atlantic postprocessing (fuel, trajectory, altitude profile).

    Only logs the columns required by the postprocessing pipeline
    (postprocessing_v1_2 / bnd_plots_v2) to minimise file size and I/O.

    Created on 2025
    @author: ahre_ni, Nils Ahrenhold
"""
import numpy as np
from bluesky import core, stack, traf, sim, navdb
from bluesky.core import Entity
from bluesky.tools import datalog, geo
from bluesky.tools.aero import ft, kts, nm, fpm
from bluesky.tools.position import txt2pos
import bluesky as bs

# ---------------------------------------------------------------------------
# Compact header — 13 columns only
# ---------------------------------------------------------------------------
flstheader = (
    'simt,'
    'flightid,'
    'ac_type,'
    'flighttime,'
    'distanceflown,'
    'totalfuel,'
    'positivefuelflow,'
    'latitude,'
    'longitude,'
    'altitude,'
    'tas,'
    'vs,'
    'spawntime\n'
)

confheader = (
    'simt[s],'
    'ac1,'
    'ac2,'
    'latitude_ac1[deg],'
    'longitude_ac1[deg],'
    'altitude_ac1[ft],'
    'latitude_ac2[deg],'
    'longitude_ac2[deg],'
    'altitude_ac2[ft],'
    'heading_ac1[deg],'
    'heading_ac2[deg],'
    'vs_ac1[fpm],'
    'vs_ac2[fpm],'
    'dcpa[nm],'
    'tcpa[sec],'
    'tLOS[sec],'
    'qdr[deg],'
    'dist[nm],'
    'dalt_min[ft],'
    'tinconf[sec],'
    'toutconf[sec],'
    'duration[s],'
    'intrusion\n'
)


# ---------------------------------------------------------------------------
# Plugin entry-point
# ---------------------------------------------------------------------------
def init_plugin():
    """Plugin initialisation function."""
    global loggercompact
    loggercompact = LoggerCompact()

    config = {
        'plugin_name': 'LOGGERCOMPACT',
        'plugin_type': 'sim',
    }

    stackfunctions = {
        'STARTLOG': [
            'STARTLOG',
            '',
            loggercompact.start_log,
            'Starts the compact flight-state and conflict logger'
        ],
        'SETMASS': [
            'SETMASS acid,mass',
            '[acid,float]',
            loggercompact.setmass,
            'Set aircraft mass [kg] after creation, e.g. SETMASS callsign, 257656.4'
        ]
    }

    return config, stackfunctions


# ---------------------------------------------------------------------------
# Logger entity
# ---------------------------------------------------------------------------
class LoggerCompact(Entity):
    """Compact logger: logs only the columns needed for postprocessing."""

    def __init__(self):
        super().__init__()

        # Destination proximity parameter [nm]
        self.d = 10

        # Conflict tracking dicts
        self.duration = {}
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
        self.dcpa = {}
        self.dalt = {}
        self.tLOS = {}
        self.dist = {}
        self.qdr = {}
        self.tcpa = {}
        self.tinconf = {}
        self.toutconf = {}

        self.sim_name = stack.get_scenname()

        # Log files
        self.flst = datalog.crelog('FLSTLOG_COMPACT', None, flstheader)
        self.conflog = datalog.crelog('CONFLOG_COMPACT', None, confheader)

        with self.settrafarrays():
            self.create_time = np.array([])
            self.total_fuel = np.array([])
            self.last_update_time = np.array([])

    # ------------------------------------------------------------------
    def reset(self):
        super().reset()
        self.duration = {}
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
        self.dcpa = {}
        self.dalt = {}
        self.tLOS = {}
        self.dist = {}
        self.qdr = {}
        self.tcpa = {}
        self.tinconf = {}
        self.toutconf = {}
        self.sim_name = None

    # ------------------------------------------------------------------
    def create(self, n=1):
        super().create(n)
        self.create_time[-n:] = sim.simt
        self.total_fuel[-n:] = 0.0
        self.last_update_time[-n:] = sim.simt

    # ------------------------------------------------------------------
    @core.timed_function(name='LOGGERCOMPACT', dt=1.0)
    def update(self, dt):
        """
        1. Handle aircraft deletion (destination proximity)
        2. Early exit if no aircraft
        3. Fuel integration
        4. Conflict processing
        5. Compact flight-state logging
        """
        # ---------------------------------------------------------------
        # 1. AIRCRAFT DELETION
        # ---------------------------------------------------------------
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
                print(f"COMPACT LOGGER - {self.sim_name}: {cs} landed at {sim.simt}; "
                      f"active aircraft: {traf.ntraf}")

            # Safety net: delete zombie aircraft that have been flying
            # longer than 24 h or have invalid coordinates.
            # This catches flights whose destination lookup failed
            # (valid_dest=False) so they were never distance-checked.
            if traf.ntraf > 0:
                flight_ages = sim.simt - self.create_time[:traf.ntraf]
                max_flight_s = 24.0 * 3600.0           # 24 hours hard limit
                bad_lat = np.abs(traf.lat[:traf.ntraf]) > 90.0
                too_old = flight_ages > max_flight_s
                zombie_mask = too_old | bad_lat
                zombie_indices = np.where(zombie_mask)[0]
                # Delete in reverse order to keep indices valid
                for idx in sorted(zombie_indices, reverse=True):
                    cs = traf.id[idx]
                    age_h = flight_ages[idx] / 3600.0
                    traf.delete(idx)
                    print(f"COMPACT LOGGER - {self.sim_name}: ZOMBIE {cs} deleted "
                          f"(age={age_h:.1f}h); active aircraft: {traf.ntraf}")

        # ---------------------------------------------------------------
        # 2. EARLY EXIT
        # ---------------------------------------------------------------
        if sim.simt >= 48 * 3600 and traf.ntraf == 0:
            name = self.sim_name or "UNNAMED"
            print(f"END of simulation: {name} at {sim.simt}s = {sim.simt/3600:.1f}h")
            stack.stack('reset')
            return

        if traf.ntraf == 0:
            return

        n = traf.ntraf

        # ---------------------------------------------------------------
        # 3. FUEL INTEGRATION
        # ---------------------------------------------------------------
        try:
            raw_ff = traf.perf.fuelflow[:n]
            pos_ff = np.maximum(0, raw_ff)
            self.total_fuel[:n] += pos_ff * dt
            self.last_update_time[:n] = sim.simt
        except (AttributeError, IndexError):
            pos_ff = np.full(n, 0.0)

        # ---------------------------------------------------------------
        # 4. CONFLICT PROCESSING
        # ---------------------------------------------------------------
        idxdict = {frozenset(v): i for i, v in enumerate(traf.cd.confpairs)}
        for cpf in traf.cd.confpairs_unique:
            if cpf in idxdict:
                i = idxdict[cpf]
                self.dcpa[cpf] = np.asarray(traf.cd.dcpa)[i]
                self.dalt[cpf] = np.asarray(traf.cd.dalt)[i]
                self.tcpa[cpf] = np.asarray(traf.cd.tcpa)[i]
                self.tLOS[cpf] = np.asarray(traf.cd.tLOS)[i]
                self.qdr[cpf] = np.asarray(traf.cd.qdr)[i]
                dist_now = np.asarray(traf.cd.dist)[i]
                if cpf not in self.dist:
                    self.dist[cpf] = dist_now
                else:
                    self.dist[cpf] = min(self.dist[cpf], dist_now)
                if cpf not in self.tinconf:
                    self.tinconf[cpf] = sim.simt
                if cpf in traf.cd.lospairs_unique:
                    self.intrusion_occurred[cpf] = True
                elif cpf not in self.intrusion_occurred:
                    self.intrusion_occurred[cpf] = False

        # Duration tracking & end-of-conflict logging
        for pair in traf.cd.confpairs:
            up = frozenset(pair)
            if up not in self.duration:
                self.duration[up] = 1
                ac1, ac2 = tuple(pair)
                i1, i2 = traf.id2idx(ac1), traf.id2idx(ac2)
                self.init_lat1[up] = traf.lat[i1]
                self.init_lon1[up] = traf.lon[i1]
                self.init_alt1[up] = traf.alt[i1]
                self.init_hdg1[up] = traf.hdg[i1]
                self.init_vs1[up] = traf.vs[i1] / fpm
                self.init_lat2[up] = traf.lat[i2]
                self.init_lon2[up] = traf.lon[i2]
                self.init_alt2[up] = traf.alt[i2]
                self.init_hdg2[up] = traf.hdg[i2]
                self.init_vs2[up] = traf.vs[i2] / fpm

        for pair, dur in list(self.duration.items()):
            cpf = frozenset(pair)
            if cpf in traf.cd.confpairs_unique:
                self.duration[pair] += 1
            else:
                ac1, ac2 = tuple(pair)
                self.conflog.log(
                    ac1, ac2,
                    self.init_lat1[cpf], self.init_lon1[cpf], self.init_alt1[cpf] / ft,
                    self.init_lat2[cpf], self.init_lon2[cpf], self.init_alt2[cpf] / ft,
                    self.init_hdg1[cpf], self.init_hdg2[cpf],
                    self.init_vs1[cpf], self.init_vs2[cpf],
                    self.dcpa[cpf] / nm, self.tcpa[cpf], self.tLOS[cpf],
                    self.qdr[cpf], self.dist[cpf] / nm, self.dalt[cpf] / ft,
                    self.tinconf[cpf], sim.simt, self.duration[pair],
                    int(self.intrusion_occurred[cpf])
                )
                del self.duration[pair]
                self.toutconf[cpf] = sim.simt

        # ---------------------------------------------------------------
        # 5. COMPACT FLIGHT-STATE LOGGING  (13 columns)
        # ---------------------------------------------------------------
        if n > 0:
            self.flst.log(
                traf.id[:n],                               # flightid
                traf.type[:n],                             # ac_type
                sim.simt - self.create_time[:n],           # flighttime [s]
                traf.distflown[:n] / nm,                   # distanceflown [nm]
                self.total_fuel[:n],                       # totalfuel [kg]
                np.maximum(0, traf.perf.fuelflow[:n]) if hasattr(traf.perf, 'fuelflow') else np.zeros(n),
                                                           # positivefuelflow [kg/s]
                traf.lat[:n],                              # latitude [deg]
                traf.lon[:n],                              # longitude [deg]
                traf.alt[:n] / ft,                         # altitude [ft]
                traf.tas[:n] / kts,                        # tas [kts]
                traf.vs[:n] / fpm,                         # vs [fpm]
                self.create_time[:n],                      # spawntime [s]
            )

    # ------------------------------------------------------------------
    def start_log(self):
        print(f"LOGGERCOMPACT: FLST and CONF logger started: {self.sim_name}")
        self.flst.start()
        self.conflog.start()
        return True, f'Compact FLST and CONF logger is ON for simulation: {self.sim_name}.'

    # ------------------------------------------------------------------
    def setmass(self, idx, mass_kg):
        """Set aircraft mass [kg] in the BADA performance model.

        Usage in scenario file:
            SETMASS callsign, mass_in_kg
        """
        if mass_kg <= 0:
            return False, f'SETMASS: Mass must be positive, got {mass_kg} kg.'

        acid = traf.id[idx]

        if not hasattr(traf.perf, 'mass') or len(traf.perf.mass) == 0:
            return False, f'SETMASS {acid}: No performance model with mass data active.'

        if hasattr(traf.perf, 'mmin') and hasattr(traf.perf, 'mmax'):
            mmin = traf.perf.mmin[idx]
            mmax = traf.perf.mmax[idx]
            if mass_kg < mmin:
                return False, f'SETMASS {acid}: Mass {mass_kg:.1f} kg below OEW ({mmin:.1f} kg).'
            if mass_kg > mmax:
                return False, f'SETMASS {acid}: Mass {mass_kg:.1f} kg exceeds MTOW ({mmax:.1f} kg).'

        old_mass = traf.perf.mass[idx]
        traf.perf.mass[idx] = mass_kg
        print(f'SETMASS {acid}: Mass set from {old_mass:.1f} kg to {mass_kg:.1f} kg')
        return True, f'SETMASS {acid}: Mass set to {mass_kg:.1f} kg (was {old_mass:.1f} kg).'
