# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 08:49:23 2016

@author: metz_is
"""



import numpy as np


''' TODO: DATA LOG'''
#import CDatalog


import os
import pandas as pd


import bluesky as bs
from bluesky import core, stack
from bluesky.core import Entity, trafficarrays
from bluesky.tools.aero import ft, kts
from bluesky.stack.cmdparser import append_commands
#from bluesky.simulation import simt

#from plugins.birds.randomize_birdies import randomize_birds
from bluesky.core.walltime import Timer
from bluesky.traffic import Traffic 
from plugins.birdtraffic import BirdTraffic

import os
dir = os.path.dirname(__file__)


def init_plugin():
    
    global bird_cdr
    
    bird_cdr = Conflict_Detection_Birds()

    
    config = {
        # The name of your plugin
        'plugin_name'      : 'BIRDCDR',
        'plugin_type'      : 'sim',
        'update'           : update,
        }
    print ("loaded birdtraf yay")
    return config




def access_plugin_module(plugin_name):
    ''' access bluesky module by name - here birdtraf
    --> can be done via core once the merge from the main was performed'''
    
    plugin_name = plugin_name.lower()
    
    plugin_module = dict(bs.core.varexplorer.varlist)[plugin_name][0]
    
    
    return plugin_module






def update():
    # do modelling here. update bird state lat,lon,alt,trk,etc
    
    # send data to bird gui (this should be last step of your things)
    #bird_traf.release_birds()
    
    ''' no manual update rate here but checking as often as feasible for
    the simulation because approaches of birds and aircraft happen very quickly
    due to high aircraft speeds. We wouldn't want to miss a collision'''
    
    # this is Isabel's one
   # bird_cdr.conflict_detection()
    
    # this is Isha's one
    bird_cdr.collision_avoidance()
    bird_cdr.log_aircraft_trajectory()

    # make sure that bs traffic elements are updated
    bs.traf.ac_collision_radius = bird_cdr.ac_collision_radius
    bs.traf.ac_collision_height = bird_cdr.ac_collision_height  
    bs.traf.ac_collision_sweep  = bird_cdr.ac_collision_sweep     
    
    



    
    return

@stack.command(name = 'LOGNAME')
def LOGNAME(filename):
    '''when we want to load a scenario'''

    # bird left the area, landed or was eaten by an aircraft
    bird_cdr.set_logname(filename) 
    
    return

@stack.command(name = 'CRE_BIRDAC')
def CRE_BIRDAC(acid, actype: str="B744", aclat: float=52., aclon: float=4., achdg: float=None, acalt: float=0,  
        acspd: float = 0):
    """CREM2 acid, type, [latlon], [hdg], [alt], [spd], prio"""
    ## DEPRECATED!!!
    # Creates an aircrft, but also assigns priority
    # Convert stuff for bs.traf.cre

    # correct some argument units
    acspd *= kts
    acalt *= ft
        
    # First create the aircraft

    bs.traf.cre(acid, actype, aclat, aclon, achdg, acalt, acspd)

       
    # Then assign its collision envelope
    idx = bs.traf.id.index(acid)

    
    coll_rad, coll_height, coll_sweep = bird_cdr.assign_envelope(actype)

    
  #  bs.traf.priority[idx] = prio

    # you can just do this
    bird_cdr.ac_collision_radius[idx] = coll_rad
    bird_cdr.ac_collision_height[idx] = coll_height
    bird_cdr.ac_collision_sweep[idx]  = coll_sweep
    
    # bs.traf.ac_collision_radius[idx] = 342.
    # bs.traf.ac_collision_height[idx] = coll_height
    # bird_cdr.traf.ac_collision_radius[idx] = 18.
    # bird_cdr.traf.ac_collision_height[idx] = coll_height


    # add path plan for specific aircraft
    #bs.traf.path_plans[-1] = path_plan_dict[ACID]

    return









##############################################################################
###
###
###
##############################################################################




class Conflict_Detection_Birds(Entity):
    
    def __init__(self):
        
        super().__init__()
        ''' get our birdie info from birdtraffic'''
        
        birdtraf_module = access_plugin_module('BIRDSIM')
        self.birds = birdtraf_module.bird_traf
        self.traf = bs.traf
        
        with self.settrafarrays():
            # initialization of extra traffic information, in this case
            #  collision envelope. Add variables for other parts of the 
            # safety envelope here
            
            self.ac_collision_radius = np.array([], dtype=int)
            self.ac_collision_height = np.array([], dtype=int)
            self.ac_collision_sweep = np.array([], dtype=int)

            # print ('we initialized', self.traf.ac_collision_radius, self.ac_collision_radius)

        # and make life with regard to traf a bit easier
        bs.traf.ac_collision_radius = self.ac_collision_radius
        bs.traf.ac_collision_height = self.ac_collision_height  
        bs.traf.ac_collision_sweep = self.ac_collision_sweep 
        


        self.counter_strikes = 0
        self.counter = 0
        
        global earth_radius
        earth_radius = 6371000.0 # Earth radius in m          
        self.coslat = np.cos(np.radians(51.4192475))
        # Create datalog instance

        self.log = Datalog()
        self.reset()
        

        
        
        
        # for logging we need a filename. If we use the IC, the name is given
        # this is a placeholder to avoid errors in case no logname is provided in the 
        # scenario
        self.filename2save = 'ZZZZ_bird_logging'
        
        
        ''' to perform CD, we need to know whether the bird enters the 
        safety envelope. For that purpose, we need to define the expansion of the 
        safety envelope. Here, the required parameters are initialized. In "add_envelope", 
        they will be assigned to newly created aircraft'''
        



        return
    
    def create(self, n=1):
        
        ''' this is to ensure that we have arrays for collision envelopes which
        we can fill with values later on '''
        super().create(n)


        self.ac_collision_radius[-n:] = 20. # unit is m
        self.ac_collision_height[-n:] = 1.42 # unit is m
        self.ac_collision_sweep[-n:]  = 24. # degrees

        
        

    def assign_envelope(self, ac_type):
        
        ### now we have hard-coded radii and heights. Add a function here
        ### to assign the values depending on the aircraft type
        # different numbers to init for checking only
        coll_rad = 42.
        coll_height = 1.49
        
        # sweep is relevant for fixed-wing aircraft to ensure that
        # only collisions to the front of the aircraft are counted
        # if you want the full circle, use 90 for your aircraft types
        coll_sweep = 26.
        
        return coll_rad, coll_height, coll_sweep

        


    
    def set_logname(self, filename):
        # assigning the actual logname for the logfile recording collision parameters.
        
        self.filename2save = filename    
        

        return
    
    

    def reset(self):

        
        self.counter = self.counter + 1

        self.filename_set = False
        
        
        return        
     
        
    def log_aircraft_trajectory(self):
            # copy logging info stuff here
        return
        
    def collision_avoidance(self)   :
        
        if (len(self.birds.id) <1 ) or (len(self.traf.id) <1):
            return        
        
        
       # print ('we are in Ishas function, yay ')
       #
        '''replace with self.birds.collision_radius and self.ac_collision_radius to ensure that you 
       have as many radii as you have traffic participants! This now is just a quick fix for a 1-1 comparison!'''
       
       
       # will be replaced with warning radius
        crit_dist_ac_in = np.array([20.])
        crit_dist_bird_in =  np.array([1.])
        
        
        # !!! homework for Isha: implement ac_warn_radius and ac_caution_radius as it is done for collision radius. Check for EVERY occurrence
        
        '''replace with these lines and change the values for collision radii where they are initialized in birdtraf.'''
        #crit_dist_ac = self.birds.collision_radius.reshape(len(self.birds.collision_radius), 1)
        #crit_dist_b = self.ac_collision_radius.reshape(1,len(self.ac_collision_radius))        

        crit_dist_ac = crit_dist_ac_in.reshape(len(crit_dist_ac_in), 1)
        crit_dist_b = crit_dist_bird_in.reshape(1,len(crit_dist_bird_in))           
        
        
        

        lat_birds, lon_birds, alt_birds, lat_aircraft, lon_aircraft, alt_aircraft = self.reshape()       
       # print ("lat bird", lat_birds, "lat aircraft", lat_aircraft)
        # first filter:lateral distance
        dxy = self.distance(np.radians(lat_birds), np.radians(lon_birds), np.radians(lat_aircraft), np.radians(lon_aircraft))
        
        
        print ('dxy: ', dxy, 'crit_dist reached ', dxy < crit_dist_ac + crit_dist_b)
        if dxy < crit_dist_ac + crit_dist_b:
            print ("ohoh")
            
        
        
        return
    
    
        
     
        
    def conflict_detection(self):

       # print "simt before CD", simt
        # only run if there are birds and aircraft
        if (len(self.birds.id) <1 ) or (len(self.traf.id) <1):
            return


        # once we have the first birdies, we would like to have a filename to store
        if not self.filename_set:
            self.birds.filename_judihui = self.filename2save
            self.filename_set = True
        # bring input to correct format
        lat_birds, lon_birds, alt_birds, lat_aircraft, lon_aircraft, alt_aircraft = self.reshape()       

        # first filter:lateral distance
        dxy = self.distance(np.radians(lat_birds), np.radians(lon_birds), np.radians(lat_aircraft), np.radians(lon_aircraft))
        #dxyqwik = self.quick_distance(np.radians(lat_birds), np.radians(lon_birds), np.radians(lat_aircraft), np.radians(lon_aircraft))
        #print simt, self.traf.id, self.birds.id, "dist",dxy, "kwik", dxyqwik, "deltaalt",abs(alt_birds - alt_aircraft), "alt_b",alt_birds, "alt_ac",alt_aircraft
        #print "lat b", lat_birds, "lon b", lon_birds, "lat ac", lat_aircraft, "lon ac", lon_aircraft
        #print "rad b", self.birds.collision_radius, "rad ac", self.traf.ac_collision_radius,  "h ac", self.traf.ac_collision_height
        #print
        
        # is this already in the dangerous area?
        # input for ac is already radius (diameter/2)
        # for birds we are fixed now: 0.5m individuals, 5m flocks

        c_rad_birds = self.birds.collision_radius.reshape(len(self.birds.collision_radius), 1)
        c_rad_ac = self.ac_collision_radius.reshape(1,len(self.ac_collision_radius))
      #  print ('dxy', dxy, 'c_rad_birds', c_rad_birds, 'c_rad_ac', c_rad_ac, 'sum', c_rad_birds + c_rad_ac)
        dangerous_dist = (dxy <= c_rad_birds + c_rad_ac)*1.
       # print ('dxy', dxy, 'c_rad_birds', c_rad_birds, 'c_rad_ac', c_rad_ac, 'sum', c_rad_birds + c_rad_ac)
        
        # only continue for bird-ac combinations where the lateral distance is too small
        if len(np.where(np.any(dangerous_dist == 1. , axis = 1) == True)[0]) > 0 and \
           len(np.where(np.any(dangerous_dist == 1. , axis = 0) == True)[0]) > 0 : 
              # print "dangerous dist"
              

               
               # filter
               alt_birds = alt_birds[np.where(np.any(dangerous_dist == 1. , axis = 1) == True)[0]]
               alt_aircraft = alt_aircraft[0][np.where(np.any(dangerous_dist == 1., axis = 0) == True)[0]]
               collision_height = self.ac_collision_height[np.where(np.any(dangerous_dist == 1. , axis = 0) == True)[0]]                
               
               # filter
               lat_birds    = lat_birds[np.where(np.any(dangerous_dist == 1. , axis = 1) == True)[0]]
               lon_birds    = lon_birds[np.where(np.any(dangerous_dist == 1. , axis = 1) == True)[0]]
               lat_aircraft = lat_aircraft[0][np.where(np.any(dangerous_dist == 1., axis = 0) == True)[0]]
               lon_aircraft = lon_aircraft[0][np.where(np.any(dangerous_dist == 1., axis = 0) == True)[0]]
               
               
               
               
               # is a list and has therefore to be converted 
               id_ac = np.array(self.traf.id)
               # used for later
               sweep = self.ac_collision_sweep[np.where(np.any(dangerous_dist == 1. , axis = 0) == True)[0]]
               hdg   = self.traf.hdg[np.where(np.any(dangerous_dist == 1. , axis = 0) == True)[0]]
               id_ac = id_ac[np.where(np.any(dangerous_dist == 1. , axis = 0) == True)[0]]
               id_bird = self.birds.id[np.where(np.any(dangerous_dist == 1. , axis = 1) == True)[0]]




               # altiutde difference: 
               # only birds in the same plane as the aircraft are interesting
               # input is already ac_height/2
               dangerous_alt = (abs(alt_birds - alt_aircraft) <= collision_height)*1.
              # print ('alt bird', alt_birds, 'alt ac',alt_aircraft,'collheight', collision_height, 'dangerous alt', dangerous_alt )
               
               # only continue if there are bird-ac combinations within 
               #dangerous distance AND in the same altitude band
    
                # only continue if any birds and aircraft are in the same altitude layer
               if len(np.where(np.any(dangerous_alt == 1. , axis = 1) == True)[0]) > 0 and \
                  len(np.where(np.any(dangerous_alt == 1. , axis = 0) == True)[0]) > 0 :
                    
                    
                    # filter
                    lat_birds    = lat_birds[np.where(np.any(dangerous_alt == 1. , axis = 1) == True)[0]]
                    lon_birds    = lon_birds[np.where(np.any(dangerous_alt == 1. , axis = 1) == True)[0]]
                    

                    lat_aircraft = lat_aircraft[np.where(np.any(dangerous_alt == 1., axis = 0) == True)[0]]
                    lon_aircraft = lon_aircraft[np.where(np.any(dangerous_alt == 1., axis = 0) == True)[0]]
                    

                    sweep        = sweep[np.where(np.any(dangerous_alt == 1. , axis = 0) == True)[0]]
                    hdg          = hdg[np.where(np.any(dangerous_alt == 1. , axis = 0) == True)[0]]
                    id_ac        = id_ac[np.where(np.any(dangerous_alt == 1. , axis = 0) == True)[0]]
                    id_bird      = id_bird[np.where(np.any(dangerous_alt == 1. , axis = 1) == True)[0]]

            
                    # bearing between bird and aircraft
                    bearing = self.bearing(np.radians(lat_aircraft), np.radians(lon_aircraft), np.radians(lat_birds), np.radians(lon_birds))
                    #print "simt", simt
                    #print "bearing", bearing, "ac pos", lat_aircraft, lon_aircraft, "bird pos", lat_birds, lon_birds
                    # top view of the aircraft: bird strikes only occurr if 
                    # they take place in the front half (end is wingtip)
                    # relative values required
                    pacman_high = ( 90. + sweep)
                    pacman_low  = (-90. - sweep)
                    #print "pacmaaaan", pacman_low, pacman_high
                    # explanation in method
                    delta_heading = ((((hdg - bearing)%360.) + 180. + 360.)% 360.) - 180.        
                   # print ("delta heading", delta_heading)
                
                    # and is it within the front area of the aircraft?
                    # then we have a strike!
                    pacman = ((delta_heading > pacman_low) & (delta_heading < pacman_high) )* 1.
                    #print "pacman", pacman
                    # which birds were hit? 

                    id_hit_birds = id_bird[np.where(np.any(pacman ==1., axis = 1) == True)[0]]
                    id_hit_ac = id_ac[np.where(np.any(pacman ==1., axis = 0) == True)[0]]
                    #print
                    #print "hdg aircraft", hdg
                    #print "pacman", pacman_high, pacman_low, "delta hdg", delta_heading, "bearing",bearing
                    
                    # only continue if there was a strike
                    if len(id_hit_birds) > 0:
                        strike_time = bs.sim.simt

                        idx_birds_hit = []
                        bird_data = []
                        lat_birds = []
                        lon_birds = []
                        
                        for identity in id_hit_birds:
                            # this is the index in the class birds
                            index_birds = int(np.where(self.birds.id == float(identity))[0][0])
                            idx_birds_hit.append(index_birds) 
                            lat_birds.append(self.birds.lat[index_birds])
                            lon_birds.append(self.birds.lon[index_birds])
                            
                            bird_data.append(str(self.birds.id[index_birds]) + ' \t ' +  str(self.birds.tas[index_birds]) \
                                              + ' \t ' + str(self.birds.lat[index_birds]) + ' \t ' +  str(self.birds.lon[index_birds]) \
                                              + ' \t ' + str(self.birds.alt[index_birds]) + ' \t ' + str(self.birds.bird_size[index_birds])\
                                              + ' \t ' + str(self.birds.collision_radius[index_birds]) + ' \t ' + str(self.birds.no_inds[index_birds]) \
                                              + ' \t ' + str(self.birds.flock_flag[index_birds]))   
                        
                            # log data
                            #self.log.write(str(strike_time), "BIRD", str(self.birds.id[index_birds]),\
                            #                str(self.birds.tas[index_birds]), str(self.birds.lat[index_birds]), \
                             #               str(self.birds.lon[index_birds]), str(self.birds.alt[index_birds]), \
                             #               str(self.birds.cat[index_birds]), str(self.birds.flock_flag[index_birds]), \
                               #             "BUFFER")
                                                                    
                            
                        
                        # remove them        
                        self.birds.remove_bird(idx_birds_hit)

          
                        # store IDs of hit aircraft
                        #id_hit_ac = id_ac[np.where(np.any(pacman ==1., axis = 0) == True)[0]]

                        
                        # increase counter
                        self.counter_strikes = self.counter_strikes + len(id_hit_ac)

                        #  store the aircraft id's of the hit aircraft - preparation
                        to_mark = []
                        for identity in id_hit_ac:

                            if identity in self.traf.id:
                                to_mark.append(int(np.where(np.array(self.traf.id) == identity)[0][0]))
                        to_mark = np.unique(to_mark)
                       
                        # store the aircraft indices of the hit aircraft - execution 
                        # idx is the index of the array: 0:n
                        # pos is the value of to_mark at the index - in this case
                        # it marks the position of the aircraft within the hit_ac array
                        
                        for idx in to_mark:

                            
                            ac_data = str(self.traf.id[idx]) + ' \t ' + str(self.traf.tas[idx]) \
                             + ' \t ' + str(self.traf.lat[idx]) + ' \t ' + str(self.traf.lon[idx]) \
                             + ' \t ' + str(self.traf.alt[idx]) + ' \t ' + str(self.traf.type[idx])
                            
                            # which bird did this aircraft hit?
                            # determination via lat-lon difference (max. 0.001 resp.)
                            for i in range(len(idx_birds_hit)):

                                if (abs(self.traf.lat[idx] - lat_birds[i]) < 0.001) and \
                                    (abs(self.traf.lon[idx] - lon_birds[i] < 0.001)):

                                        #self.log.write(self.birds.filename2save, str(strike_time), ac_data, bird_data[i], "collision")
                                        self.log.write(self.birds.filename_judihui, str(strike_time), ac_data, bird_data[i], "collision")


                                        
                            # log data
                           # self.log.write(str(strike_time), "AIRCRAFT", str(self.traf.id[idx]), str(self.traf.tas[idx]), \
                            #               str(self.traf.lat[idx]), str(self.traf.lon[idx]), str(self.traf.alt[idx]), str(self.traf.type[idx]), \
                             #              str(self.traf.orig[idx]), str(self.traf.dest[idx]))
                                           
                         # write: str(strike_time), aircraft, bird --> header true/false?



                        # store the aircraft id's as well
                        # even if an aircraft has more than one strike: we only want to store it's id once
                        # np.in1d: is arg1 in arg2?
                        #new_id_hit_ac = id_hit_ac[np.where(np.in1d( id_hit_ac, self.traf.nr_strikes, invert = True))]        
                        #self.traf.nr_strikes = np.append(self.traf.nr_strikes, new_id_hit_ac)
                        #print "traf_strikes", self.traf.nr_strikes

                    # log data of collision
                        

                        #self.log.save(self.birds.filename2save) 
                        self.log.save(self.birds.filename_judihui)
                        print('saved', self.birds.filename_judihui)





        
    
        return 



    # format input for calculation
    # hint: name the reshapes differently oooooooor make an individual module
    # individual module might have adavantages as there are inputs from traf. and from birds.
    # height and radius of aircraft: store in input_files or find ways to get it via other parameters
    def reshape(self):


        # birds are the columns, aircraft are the rows
        lat_birds = self.birds.lat.reshape((len(self.birds.lat),1))
        lon_birds = self.birds.lon.reshape((len(self.birds.lon),1))
        alt_birds = self.birds.alt.reshape((len(self.birds.alt),1))

        
        lat_aircraft = self.traf.lat.reshape((1,len(self.traf.lat)))
        lon_aircraft = self.traf.lon.reshape((1,len(self.traf.lon)))
        alt_aircraft = self.traf.alt.reshape((1,len(self.traf.alt)))        
        
        
        
        return lat_birds, lon_birds, alt_birds, lat_aircraft, lon_aircraft, alt_aircraft
        
        
    # use the haversine formula to calculate the distance between birds and ac    
    # input is already in radians
    def distance(self, lat_birds, lon_birds, lat_ac, lon_ac):


        
        a = np.sin((lat_birds-lat_ac)/2)*np.sin((lat_birds-lat_ac)/2) + \
            np.cos(lat_birds)*np.cos(lat_ac)*np.sin((lon_birds-lon_ac)/2)*np.sin((lon_birds-lon_ac)/2)
   
        c= 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))
        distance= earth_radius*c # 6317000m corresponds to the earth radius

        
        return distance


    def quick_distance(self,lat_birds, lon_birds, lat_ac, lon_ac):
        '''a bit less accurracy but sooooo much faster '''
        dx = earth_radius * (lon_birds - lon_ac) * self.coslat
        dy = earth_radius * (lat_birds - lat_ac)
        distance = np.sqrt(dx*dx + dy*dy)
        
        return distance



    def bearing(self, lat1, lon1, lat2, lon2):
    
        deltal = lon2-lon1
    
    # calculate runway bearing
        bearing = np.arctan2(np.sin(deltal)*np.cos(lat2), (np.cos(lat1)*np.sin(lat2)-
                np.sin(lat1)*np.cos(lat2)*np.cos(deltal)))
        
        # normalize to 0-360 degrees
        bearing = (np.degrees(bearing)+360)%360
        
        return bearing


class Datalog():
    def __init__(self):
        print ("we are in CDR datalog")
# Create a buffer and save filename

        self.buffer=[]
        
        #filename will be set in first run
        self.filename_flag = False

         
        return
    
    def write(self, filename,time, ac_data, bird_data, occurrence_type):


        # filename[5:15] is the date
        self.buffer.append( filename +" \t "  + time +" \t " + ac_data + " \t " + bird_data + '\t' + occurrence_type + chr(13) + chr(10))
       
        return

    def save(self, filename):
        
        # files are saved per airport. Hence only create a new file if 
        # no file for this airport exists
        log_file = os.path.join(dir, "bird_CDR_log/"  + filename + ".txt" )
        print ('in save, logpath is', log_file)
        if not os.path.isfile(log_file):  
           # log_file = "log/" + filename_def + ".txt"
           # self.log_file = os.path.join(dir, log_file)
            #print "INIT", filename_def, filename, log_file
            
            with open(log_file, "a") as writeto:
                writeto.write('date \t time \t id_ac \t tas \t lat \t lon \t alt \t type \t id_bird \t tas \t lat \t lon \t alt \t size \t coll_rad \t number \t flock_flag \t occurrence type \n')
        
        
       # if not self.filename_flag:
          #  filename = "log/" + filename + ".txt"
            
           # self.filename = os.path.join(dir,  filename)

            # write the header            
            #with open(self.filename, "a") as writeto:
            #    writeto.write('time \t id_ac \t tas \t lat \t lon \t alt \t type \t orig \t dest id_bird \t tas \t lat \t lon \t alt \t size \t coll_rad \t number \t flock_flag \n')

            #self.filename_flag = True

# Write buffer to file 

        with open(log_file, "a") as writeto:
            for i in range(len(self.buffer)):

                writeto.write(self.buffer[i])


        self.buffer = []    
        


        
        
        
        
        return
   
