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
from bluesky.tools.aero import ft,kts,nm,fpm
from bluesky.tools.position import txt2pos
import bluesky as bs
from bluesky.traffic.performance.perfbase import PerfBase

flstheader = \
    'simt,' + \
    'callsign,' + \
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
    'vs,' + \
    'heading,' + \
    'asasactive,' + \
    'pilotalt,' + \
    'pilottas,' + \
    'pilothdg,' + \
    'pilotvs,' + \
    'n_active_conflicts,' + \
    'n_active_intrusions' + '\n'

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

        self.mass_warning_callsigns = set()  # Track which aircraft IDs we've warned about

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

        self.mass_warning_callsigns = set()  # Reset mass warning tracking

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
            deleted_callsign = traf.id[idx]
            self.mass_warning_callsigns.discard(deleted_callsign)
            
            # Optional: Log cleanup for debugging
            # print(f"LOGGERFF: Cleaned up tracking for deleted aircraft: {deleted_callsign}")

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
        if sim.simt >= 24 * 3600:  # 24 hours in seconds
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
        ################# 3. UPDATE TRACKING ARRAYS ##########
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

            # Thrust
            thrust = np.zeros(n_current)
            if hasattr(traf.perf, 'thrust') and len(traf.perf.thrust) >= n_current:
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

        # Ensure ALL arrays have exactly n_current elements and log
        if n_current > 0:
            self.flst.log(
                traf.id[:n_current],                               # [n_current]
                traf.type[:n_current],                             # [n_current]
                (sim.simt - self.create_time[:n_current]),
                current_mass,                                       # [n_current]
                (traf.distflown[:n_current] / nm),                 # [n_current]
                self.total_fuel[:n_current],                       # [n_current]
                traf.lat[:n_current],                              # [n_current]
                traf.lon[:n_current],                              # [n_current]
                self.create_time[:n_current],                      # [n_current]
                self.distance2D[:n_current],                       # [n_current]
                self.distance3D[:n_current],                       # [n_current]
                (traf.work[:n_current] * 1e-6),                    # [n_current]
                positive_fuelflow,                                 # [n_current]
                raw_fuelflow,
                thrust,                                            # [n_current]
                (traf.alt[:n_current] / ft),                       # [n_current]
                (traf.tas[:n_current] / kts),                      # [n_current]
                (traf.vs[:n_current] / fpm),                       # [n_current]
                traf.hdg[:n_current],                              # [n_current]                          # [n_current] - Destination Lon
                traf.cr.active[:n_current],                        # [n_current]
                (traf.aporasas.alt[:n_current] / ft),              # [n_current]
                (traf.aporasas.tas[:n_current] / kts),             # [n_current]
                traf.aporasas.hdg[:n_current],                     # [n_current]
                (traf.aporasas.vs[:n_current] / fpm),              # [n_current]
                n_active_conflicts,                                # Scalar
                n_active_intrusions                                # Scalar
            )
            
            
    def start_log(self):
        print (f"LOGGERFF: FLST and CONF loggerff started: {self.sim_name}")
        self.flst.start()
        self.conflog.start()
        return True, f'FLST and CONF loggerff is ON for simulation: {self.sim_name}.'
