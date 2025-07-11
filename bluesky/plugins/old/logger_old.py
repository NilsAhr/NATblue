""" BlueSky logger_old for conflict and flight state in North Atlantic
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
    'pilotvs' + '\n'

# header for total count of conflict
confheader = \
    'simt[s],' + \
    'totalnumberofconflicts\n'

# header for conflict parameter
#dcpa = Distance at closest point of approach [nm]
#tcpa = Time at closest point of approach [sec]
#tLOS = Time until horizontal loss of separation [sec]
#qdr  =
#dist = 
confparamheader = \
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
    'tinconf[sec],' + \
    'toutconf[sec],' + \
    'dalt_min[ft]\n'

confseverheader = \
    'simt[s],' + \
    'ac1,' + \
    'ac2,' + \
    'duration[s]\n'


### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
def init_plugin():
    ''' Plugin initialisation function. '''
    # Instantiate our example entity
    global logger_old
    logger_old = Logger_old()

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'LOGGER_old',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim',
        }
    
    stackfunctions = {
        'STARTLOG': [
            'STARTLOG',
            '',
            logger_old.start_log,
            'Starts the flight status and conflict logger_old'
        ]
        }

    # init_plugin() should always return a configuration dict.
    return config, stackfunctions


### Entities in BlueSky are objects that are created only once (called singleton)
### which implement some traffic or other simulation functionality.
### To define an entity that ADDS functionality to BlueSky, create a class that
### inherits from bluesky.core.Entity.
### To replace existing functionality in BlueSky, inherit from the class that
### provides the original implementation (see for example the asas/eby plugin).
class Logger_old(Entity):
    ''' Example new entity object for BlueSky. '''
    def __init__(self):
        super().__init__()
        # All classes deriving from Entity can register lists and numpy arrays
        # that hold per-aircraft data. This way, their size is automatically
        # updated when aircraft are created or deleted in the simulation.

        # Parameters for conflict count and statistics
        self.prevconfpairs = set()      #counter for all conflicts already listed
        self.conf_all = 0               #counter for all conflicts
        self.duration = {}              #dict for duration values
        self.d = 10                     #distance parameter for deleting aircraft automatically 

        # new severity parameters
        self.min_dcpa = {}    # Track minimum dcpa for each conflict
        self.min_dalt = {}    # Track minimum dalt for each conflict
        self.tinconf = {}     # Track conflict start times for each conflict
        self.toutconf = {}    # Track conflict end times for each conflict
        self.sim_name = stack.get_scenname()

        # The FLST & CONF & CONFPARAM logger_old
        self.flst = datalog.crelog('FLSTLOG_LOGGER', None, flstheader)
        self.conflog = datalog.crelog('CONFLOG_LOGGER', None, confheader)
        self.confparamlog = datalog.crelog('CONFPARAMLOG_LOGGER', None, confparamheader)
        self.confseverlog = datalog.crelog('CONFSEVERLOG_LOGGER', None, confseverheader)

        with self.settrafarrays():
            self.distance2D = np.array([])
            self.distance3D = np.array([])
            self.create_time = np.array([])

    def reset(self):
        super().reset()
        self.prevconfpairs = set()
        self.conf_all = 0
        self.duration = {}

        # severity parameters
        self.min_dcpa = {}
        self.min_dalt = {}
        self.tinconf = {}
        self.toutconf = {}
        self.sim_name = None


    def create(self, n=1):
        ''' This function gets called automatically when new aircraft are created. '''
        # Don't forget to call the base class create when you reimplement this function!
        super().create(n)        
        self.create_time[-n:] = sim.simt

    # Functions that need to be called periodically can be indicated to BlueSky
    # with the timed_function decorator
    @core.timed_function(name='LOGGER_old', dt=1.0)
    def update(self, dt):
        ''' Periodic update function for our example entity. '''
         
        resultantspd = np.sqrt(traf.gs * traf.gs + traf.vs * traf.vs)
        self.distance2D += dt * traf.gs
        self.distance3D += dt * resultantspd
        
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
                


        # vectorized version
        # Extract the latitudes and longitudes of all aircraft
        csn = np.array(traf.id)           # All aircraft callsign
        lat_ac = np.array(traf.lat)       # All aircraft latitudes
        lon_ac = np.array(traf.lon)       # All aircraft longitudes
        dest_names = np.array(traf.ap.dest)  # Destination airport names for each aircraft

        # Initialize destination arrays
        lat_dest = np.zeros(traf.ntraf)
        lon_dest = np.zeros(traf.ntraf)

        # Boolean array to track which aircraft have valid destinations
        valid_dest = np.ones(traf.ntraf, dtype=bool)  # Start by assuming all destinations are valid
        
        
        #########################################################
        ################## DELETE IF LANDED #####################
        #########################################################
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
        #check for new conflict pairs
        confpairs_new = list(set(traf.cd.confpairs) - self.prevconfpairs)
        if confpairs_new:
            processed_conflicts = set()
            idxdict = {v: i for i, v in enumerate(traf.cd.confpairs)}

            for pair in confpairs_new:
                confpair_frozen = frozenset(pair)
                if confpair_frozen not in processed_conflicts:
                    i = idxdict[pair]
                    
                    dcpa_new = np.asarray(traf.cd.dcpa)[i]
                    tcpa_new = np.asarray(traf.cd.tcpa)[i]
                    tLOS_new = np.asarray(traf.cd.tLOS)[i]
                    qdr_new = np.asarray(traf.cd.qdr)[i]
                    dist_new = np.asarray(traf.cd.dist)[i]
                    dalt_new = abs(traf.alt[traf.id2idx(pair[0])] - traf.alt[traf.id2idx(pair[1])])

                    ac1, ac2 = pair
                    idx1 = traf.id2idx(ac1)
                    idx2 = traf.id2idx(ac2)
                    
                    lat_1, lon_1, alt_1, hdg_1, vs_1 = traf.lat[idx1], traf.lon[idx1], traf.alt[idx1], traf.hdg[idx1], traf.vs[idx1]/fpm
                    lat_2, lon_2, alt_2, hdg_2, vs_2 = traf.lat[idx2], traf.lon[idx2], traf.alt[idx2], traf.hdg[idx2], traf.vs[idx2]/fpm
                    
                    # Initialize tracking values if the conflict is new
                    if confpair_frozen not in self.min_dcpa:
                        self.min_dcpa[confpair_frozen] = dcpa_new
                        self.min_dalt[confpair_frozen] = dalt_new
                        self.tinconf[confpair_frozen] = sim.simt
                    else:
                        # Update minimum dcpa and dalt if the current values are lower
                        self.min_dcpa[confpair_frozen] = min(self.min_dcpa[confpair_frozen], dcpa_new)
                        self.min_dalt[confpair_frozen] = min(self.min_dalt[confpair_frozen], dalt_new)

                    # Log the parameters
                    self.confparamlog.log(
                        ac1, ac2, lat_1, lon_1, alt_1/ft, lat_2, lon_2, alt_2/ft, hdg_1, hdg_2, vs_1, vs_2,
                        self.min_dcpa[confpair_frozen] / nm, tcpa_new, tLOS_new, qdr_new, dist_new / nm,
                        self.tinconf[confpair_frozen], sim.simt, self.min_dalt[confpair_frozen] / ft
                    )
                    
                    processed_conflicts.add(confpair_frozen)
                    self.conf_all += 1
                    self.conflog.log(self.conf_all)

        self.prevconfpairs = set(traf.cd.confpairs)
        
        '''
        confpairs_new = list(set(traf.cd.confpairs) - self.prevconfpairs)
        if confpairs_new:

            # Dictionary to track unique conflict pairs
            processed_conflicts = set()
    
            # Mapping conflict pairs to their index in the original list
            idxdict = {v: i for i, v in enumerate(traf.cd.confpairs)}
    
            
            for pair in confpairs_new:
                # Use frozenset to ensure the conflict pair is considered in a unique way (order-independent)
                confpair_frozen = frozenset(pair)
        
                if confpair_frozen not in processed_conflicts:
                    # Get the index of the conflict pair
                    i = idxdict[pair]
            
                    # Retrieve parameters for the conflict pair
                    dcpa_new = np.asarray(traf.cd.dcpa)[i]
                    tcpa_new = np.asarray(traf.cd.tcpa)[i] #Time to closest point of approach (CPA) between aircraft
                    tLOS_new = np.asarray(traf.cd.tLOS)[i]
                    qdr_new = np.asarray(traf.cd.qdr)[i]
                    dist_new = np.asarray(traf.cd.dist)[i]
                    
                    # Aircraft IDs
                    ac1, ac2 = pair
                    
                    # Get the indices for the aircraft involved in the conflict
                    idx1 = traf.id2idx(ac1)
                    idx2 = traf.id2idx(ac2)
                    
                    # Get positions for both aircraft
                    lat_1, lon_1, alt_1, hdg_1, vs_1 = traf.lat[idx1], traf.lon[idx1], traf.alt[idx1], traf.hdg[idx1], traf.vs[idx1]/fpm
                    lat_2, lon_2, alt_2, hdg_2, vs_2 = traf.lat[idx2], traf.lon[idx2], traf.alt[idx2], traf.hdg[idx2], traf.vs[idx2]/fpm
                    
                    # Log conflict parameters
                    self.confparamlog.log(
                        ac1, ac2, lat_1, lon_1, alt_1/ft, lat_2, lon_2, alt_2/ft, hdg_1, hdg_2, vs_1, vs_2,
                        dcpa_new / nm, tcpa_new, tLOS_new, qdr_new, dist_new / nm
                    )
            
                    # Mark the conflict as processed
                    processed_conflicts.add(confpair_frozen)
                    
                    # Optional: Count and log all unique conflicts
                    self.conf_all += 1
                    self.conflog.log(self.conf_all)

        # Update the set of previous conflict pairs
        self.prevconfpairs = set(traf.cd.confpairs)
        '''

        ##########################################################
        ################## CONFLICT DURATION #####################
        ########################################################## 
        # 
        # - conflicts still ongoing when simulation ends are not logged
        #       
        # Loop through each conflict pair in traf.cd.confpairs
        for pair in traf.cd.confpairs:
        # If the conflict pair is new, initialize its duration to 1
            # Use frozenset to store unique pairs in duration dictionary
            unique_pair = frozenset(pair)
            if unique_pair not in self.duration:
                self.duration[unique_pair] = 1
            else:
                # If the conflict pair already exists, its duration has already started tracking
                continue

        # Now, loop through each pair in the duration dictionary
        for pair, duration in list(self.duration.items()):  # Use list() to allow modifying the dict during iteration
            confpair_frozen = frozenset(pair)

            # Check if the pair is still part of traf.cd.confpairs_unique
            if confpair_frozen in traf.cd.confpairs_unique:
                # Conflict is still ongoing, increment the duration by 1
                self.duration[pair] += 1
            else:
                # Conflict has ended, log the pair and its duration
                ac1, ac2 = tuple(pair) # Unpack frozenset to individual aircraft IDs
                self.confseverlog.log(ac1, ac2, self.duration[pair])

                # Remove the pair from the duration dictionary as it is no longer active
                del self.duration[pair]
                self.toutconf[confpair_frozen] = sim.simt  # Log conflict end time

        

        ################## FLIGHT STATUS LOG #####################
        #log flight statistics
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
            )
            
            

    def start_log(self):
        print (f"LOGGER: FLST and CONF logger_old started: {self.sim_name}")
        self.flst.start()
        self.conflog.start()
        self.confparamlog.start()
        self.confseverlog.start()
        return True, f'FLST and CONF logger_old is ON for simulation: {self.sim_name}.'


    # You can create new stack commands with the stack.command decorator.
    # By default, the stack command name is set to the function name.
    # The default argument type is a case-sensitive word. You can indicate different
    # types using argument annotations. This is done in the below function:
    # - The acid argument is a BlueSky-specific argument with type 'acid'.
    #       This converts callsign to the corresponding index in the traffic arrays.
    # - The count argument is a regular int.
    #@stack.command
    #def passengers(self, acid: 'acid', count: int = -1):
     #   ''' Set the number of passengers on aircraft 'acid' to 'count'. '''
      #  if count < 0:
       #     return True, f'Aircraft {traf.id[acid]} currently has {self.npassengers[acid]} passengers on board.'

        #self.npassengers[acid] = count
        #return True, f'The number of passengers on board {traf.id[acid]} is set to {count}.'
