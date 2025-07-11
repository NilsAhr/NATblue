""" BlueSky logger for conflict and flight state in North Atlantic
Created on 2024 May
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

flstheader = \
    'simt,' + \
    'callsign,' + \
    'ac_type,' + \
    'spawntime,' + \
    'flighttime,' + \
    'distanceflown,' + \
    'actualdistance2D,' + \
    'actualdistance3D,' + \
    'workdone,' + \
    'latitude,' + \
    'longitude,' + \
    'altitude,' + \
    'tas,' + \
    'vs,' + \
    'heading,' + \
    'originlat,' + \
    'originlon,' + \
    'destinationlat,' + \
    'destinationlon,' + \
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
    global logger
    logger = Logger()

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'LOGGER',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim',
        }
    
    stackfunctions = {
        'STARTLOG': [
            'STARTLOG',
            '',
            logger.start_log,
            'Starts the flight status and conflict logger'
        ]
        }

    # init_plugin() should always return a configuration dict.
    return config, stackfunctions


class Logger(Entity):
    ''' Example new entity object for BlueSky. '''
    def __init__(self):
        super().__init__()

        # Parameters for conflict count and statistics
        self.duration = {}              #dict for duration values
        self.d = 10                     #distance parameter for deleting aircraft automatically 

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

        # The FLST & CONF LOGGER
        self.flst = datalog.crelog('FLSTLOG_LOGGER', None, flstheader)
        self.conflog = datalog.crelog('CONFLOG_LOGGER', None, confheader)

        with self.settrafarrays():
            self.distance2D = np.array([])
            self.distance3D = np.array([])
            self.create_time = np.array([])

    def reset(self):
        super().reset()
        self.duration = {}

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
        # Don't forget to call the base class create when you reimplement this function!
        super().create(n)        
        self.create_time[-n:] = sim.simt

    @core.timed_function(name='LOGGER', dt=1.0)
    def update(self, dt):
         
        resultantspd = np.sqrt(traf.gs * traf.gs + traf.vs * traf.vs)
        self.distance2D += dt * traf.gs
        self.distance3D += dt * resultantspd
        
        #################################################
        ##### STOP SIMULATION IF NO ACTIVE AIRCRAFT #####
        #################################################

        # Check after 24 hours of simulation time if there are active aircraft
        if sim.simt >= 24 * 3600:  # 24 hours in seconds
            if traf.ntraf == 0:
                if self.sim_name is None:
                    print(f"END of simulation: NAME NOT DEFINED at: {sim.simt}seconds = {sim.simt/3600}hours")
                else:
                    print(f"END of simulation: {self.sim_name} at: {sim.simt}seconds = {sim.simt/3600}hours")
                stack.stack('hold')
                # stop would terminate the simulation immediately and batch simulations would stop too
                #stack.stack('stop')
                

        #########################################################
        ################## DELETE IF LANDED #####################
        #########################################################

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
            print(f"FLST LOGGER - {self.sim_name}: {cs} landed at {sim.simt}; active aircraft: {ac_act}")
        
        
        #################################################################
        ################## CONFLICT PAIRS PARAMETER #####################
        #################################################################        
        
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

        
        ##########################################################
        ################## LOG FLIGHT STATUS #####################
        ##########################################################
        # Count active conflicts and intrusions
        n_active_conflicts = len(traf.cd.confpairs_unique)
        n_active_intrusions = len(traf.cd.lospairs_unique)

        self.flst.log(
            traf.id, #callsign
            traf.type, #aircraft type
            self.create_time, #creation time [s]
            sim.simt - self.create_time, #flight time [s]
            traf.distflown/nm, #distance flown [nm]
            self.distance2D, #actual distance 2D [nm]
            self.distance3D, #actual distance 3D [nm]
            traf.work*1e-6, #work done [MJ]
            traf.lat, #latitude [deg]
            traf.lon, #longitude [deg]
            traf.alt/ft, #altitude [ft]
            traf.tas/kts, #TAS [kts]
            traf.vs/fpm, #Vertical speed [fpm]
            traf.hdg, #heading [deg]
            "",#navdb.aptlat[traf.ap.orig], #Origin Lat [deg]
            "",#navdb.aptlon[traf.ap.orig], #Origin Lon [deg]
            "",#navdb.aptlat[traf.ap.dest], #Destination Lat [deg]
            "",#navdb.aptlon[traf.ap.dest], #Destination Lon [deg]
            traf.cr.active, #ASAS active [bool]
            traf.aporasas.alt/ft, #Pilot alt [ft]
            traf.aporasas.tas/kts, #Pilot SPD (TAS) [kts]
            traf.aporasas.hdg, #Pilot HDG [deg]
            traf.aporasas.vs/fpm, #Pilot vertical speed [fpm]
            n_active_conflicts, # Number of active conflicts
            n_active_intrusions # Number of active intrusions
        )
            
            
    def start_log(self):
        print (f"LOGGER: FLST and CONF logger started: {self.sim_name}")
        self.flst.start()
        self.conflog.start()
        return True, f'FLST and CONF logger is ON for simulation: {self.sim_name}.'
