''' Bird traffic simulation plugin '''
import numpy as np
import os
import pandas as pd


import bluesky as bs
from bluesky import core, stack
from bluesky.tools.aero import ft, kts
from bluesky.stack.cmdparser import append_commands
#from bluesky.simulation import simt

from plugins.randomize_birdies import randomize_birds
from bluesky.core.walltime import Timer


def init_plugin():
    
    global bird_traf
    # initialize bird traffic
    bird_traf = BirdTraffic()

    
    config = {
        # The name of your plugin
        'plugin_name'      : 'BIRDSIM',
        'plugin_type'      : 'sim',
        'update'           : update,
        'reset'            : reset
        }
    print ("loaded birdtraf yay")
    return config

def update():
    ''' not necessary anymore due to timer in init of BirdTraffic'''
    # do modelling here. update bird state lat,lon,alt,trk,etc
    
    # send data to bird gui (this should be last step of your things)
    #bird_traf.release_birds()
    return
    
    

def reset():
  #  print ('BT - call reset')
    # clear everything. TODO: smarter way to do this
    #print('STACK - call bird reset')
    bird_traf.reset()
    
    # release birds with no info to clear screen TODO: smarter way to do this
    #bird_traf.release_birds()

@stack.command(name='CREBIRD')
def CREBIRD(birdid, bird_object: str="individual", bird_size: int = 4, no_inds: int = 1, birdlat: float=52., birdlon: float=4., birdtrk: float=None, birdalt: float=0,  
        birdspd: float = 0):
    '''for individual birdies'''
    
    ''' CREBIRD birdid,fl_of_ind, size, number,lat,lon,hdg,alt,spd '''
    # correct some argument units
    if bird_object == 'individual':
        flock_flag = False
    else: flock_flag = True
    
    
    birdspd *= kts
    birdalt *= ft

    # create the bird
    bird_traf.create_individual(birdid, flock_flag, bird_size, no_inds, birdlat, birdlon, birdtrk, birdalt, birdspd)



@stack.command(name = 'DELBIRD')
def DELBIRD(birdid):
    # bird left the area, landed or was eaten by an aircraft
    
    # remove_bird needs an array index - convert
    index_to_remove = bird_traf.id2idx(birdid)
    #print ('remove from stack')
    bird_traf.remove_bird(index_to_remove)
    
@stack.command(name = 'BIRDS')
def BIRDS(filename):
    '''when we want to load a scenario'''
    # bird left the area, landed or was eaten by an aircraft
    bird_traf.crescenario(filename)    


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################      


# CLASS BIRD TRAFFIC


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################


class BirdTraffic():
   
    def __init__(self):

        
        # to find bird files, we need to know the current directory
        self.dir = os.path.dirname(__file__)
        
        # go back to the roots
      #  print('BT - call from init')
        self.reset()
        self.scenario_loaded = False
        self.id = []
        
        
        # some global parameters
        global earth_radius
        earth_radius = 6371000.0 # Earth radius in m  
        
        global individual_radius
        global flock_radius
        
        individual_radius = 0.5
        flock_radius = 5.
        
        global lat_south
        global lat_north
        global lon_west
        global lon_east
        

        lat_south = 51.2855167
        lat_north = 51.552975

        lon_west = 5.1347806
        lon_east = 5.50
        
        
         
                 # see screenio.py
        # Update rate of aircraft update messages [Hz]
        ''' bird update is called from here at predefined rate'''
        birdupdate_rate = 1

        # create a timer to send bird data
        self.fast_timer = Timer()
        self.fast_timer.timeout.connect(self.release_birds)
        self.fast_timer.start(int(1000 / birdupdate_rate))
        
        bs.traf.birdid = self.id



######################################################################################################
######################################################################################################      


# CREATE BIRDS FROM SCENARIO FILE


######################################################################################################
######################################################################################################


      
    def crescenario(self,filename):
        print ("whooooooooooooooooooooooop ", filename)
        
        self.scenario_loaded = True
        
    
        #print "in readin", datetime.datetime.now() - self.starttime, datetime.datetime.now()
    
        # required for Monte-Carlo-Sims: the full name of the file (e.g. EGKK_2016_06_03-1)
        # is needed for the logging. But as we have the same file as input for all
        # of the scenarios, only EGKK_2016_06_03 should be used for the path
        # where to read the file from
        self.filename2save = filename # used for recording
    
        
    
        
        # we do not always do MC simulations, so here is the back door
        # for the aircraft-MC-simulations we need to make sure not to run birdie MC sims
        # format of these files is EGKK_1-100-100-high/medium/low
        if (len(self.filename2save) > 15) and (self.filename2save[-1].isdigit()):
            print ("MONTE CARLO SIMULATION!", filename, filename[0:15])
            filename2use = "bird_movements/" + filename[0:15] + ".csv"
            filename_path = os.path.join(self.dir, filename2use)
            print (filename_path)
            
            # in case of MC simulations, we want randomized speed, heading and position
    
            # if the filename works: pass on to randomizer
            if os.path.isfile(filename_path):
                
                
                # seed depends on bird movement plan number
                # filename2use[15] is the underline. Hence the number starts at position filename2use[16]
                seed = int(filename[16:])
                print ("bird movement plan exists for that day, seed is ", seed)
                data = randomize_birds(seed, filename_path)
    
                self.assign_values(data)
                
                
            else:
                print ("no such MC birdie file ", filename )     
            
            
            
            
        else:

            filename2use = "bird_movements/" + filename + ".csv"
            filename = os.path.join(self.dir, filename2use)
    
            print ("we have a birdie file, we read from", filename)
        
            # if not: tell the user
            #try:
                
                    # cat means bird size
        
                    # because of the MC simulations, more  columns are required. However, we only need the limited set to continue
            data = pd.read_csv(filename, sep = "\,", \
                                 names = ['id', 'date', 'lon','lat', 'alt', 'cat', 'no_individuals', 'flock_flag',\
                                 'id1', 'hdg', 'spd', 'lat_s1', 'lat_n1', 'lon_w1', 'lon_e1','lat_s2', 'lat_n2', 'lon_w2', 'lon_e2'], index_col = False, engine = 'python')
            
                
 
            '''HAS TO GO AWAY AGAIN'''
            #data = pd.read_csv(filename, sep = "\,", \
             #    names = ['id', 'date', 'lon','lat', 'alt', 'cat',  'flock_flag',\
              #   'id1', 'hdg', 'spd', 'lat_s1', 'lat_n1', 'lon_w1', 'lon_e1','lat_s2', 'lat_n2', 'lon_w2', 'lon_e2'], index_col = False, engine = 'python')
            #data['no_individuals'] = 1
            
            '''HAS TO GO AWAY AGAIN'''

            data = data.drop(['lat_s1', 'lat_n1', 'lon_w1', 'lon_e1','lat_s2', 'lat_n2', 'lon_w2', 'lon_e2'], axis=1)
            
            # for the "next_ts", lat, lon, alt, we need to shift the respective cells
            # BUT FIRST SORT CORRECTLY!
            data = data.sort_values(['id','date'], ascending=[True, True])

            data['timeshift'] = data['date'].shift(-1)
            data['latshift'] = data['lat'].shift(-1)
            data['lonshift'] = data['lon'].shift(-1)

            
            
            '''to prevent birdies from being deleted one timestep too early respective start flying in a weird direction,
            we give them a theoretical final position after which they are deleted '''
            #print "init"
            #print data

            # the NaNs are the fly-outs - give them an extra 10 minutes
            data['timeshift'].fillna(data['date'] + 600., inplace = True)

            



            #print 
            #print "after timeshift filled"
            #print data

            # we need delta time, delta s
            data['delta_t'] = data['timeshift'] - data['date']
            
            
            data['delta_s'] = data['delta_t'] * data['spd']

            #print 
            #print "after delta s and delta t"
            #print data
            
            # and now calculate the lat and lon out
           # print ('we send to leo')
            data['lat_leo'], data['lon_leo'] = self.calculate_leo_position(data['delta_t'], data['delta_s'], data['lat'], data['lon'], data['hdg'])
            #print (data)
            #print 
            #print "after latlon_Leo"
            #print data
             
            # and now replace all nans or 0s of the df in the latshift and lonshift 
            data.loc[data['latshift'].isnull(),'latshift'] = data['lat_leo']
            data.loc[data['lonshift'].isnull(),'lonshift'] = data['lon_leo']
            
            data.loc[data['latshift']== 0,'latshift'] = data['lat_leo']
            data.loc[data['lonshift'] == 0,'lonshift'] = data['lon_leo']            
            
            #print 
            #print "after all"
            #print data
            data = data.sort_values(by='date')
            
            #why not directly convert?
            #--> INPUT IS IN DEGREES AND DUE TO NEW UPDATING WE DO NOT NEED
                # RADIANS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               


                
                
            self.assign_values(data)
                    
           # except:
            #    print ("no such individual birdie file")
            #    return
    
        
        return


        

######################################################################################################
######################################################################################################      


# UPDATE SCENARIO VALUES


######################################################################################################
######################################################################################################    

    def update_scen(self):

        # only process when there is at least one bird left

        
        if len(self.input_time) <1 or bs.sim.simt > self.input_time[-1]:
            '''do we really want to reset here? or just return?'''
            self.reset()
            return
        
        # work with all values corresponding to timestamps already over
        '''
        This is a little faster (and looks nicer)

        np.argmax(aa>5)

        Since argmax will stop at the first True ("In case of multiple occurrences of the 
        maximum values, the indices corresponding to the first occurrence are returned.") 
        and doesn't save another list.
        --> use  index in sense of input_time = self.input_time[0:index]
        '''    
            
            
        #idx_time_passed = np.where(self.input_time <= simtime)[0]  
        idx_time_passed = np.argmax(self.input_time > bs.sim.simt)
       # print "idx time passed", idx_time_passed      
        
            

        

        # bird info to the check-function
        # check bird is performed for all the birds that we got in the lists
        #print "we check", input_id1
        if idx_time_passed > 0:
            input_time            = self.input_time[ : idx_time_passed]
            input_id1             = self.input_id1[ : idx_time_passed]
            input_id2             = self.input_id2[ : idx_time_passed]
            input_bird_size       = self.input_bird_size[ : idx_time_passed]
            input_no_inds         = self.input_no_inds[ : idx_time_passed]
            input_flock_flag      = self.input_flock_flag[ : idx_time_passed]
            input_alt             = self.input_alt[ : idx_time_passed]      
            input_lat             = self.input_lat[ : idx_time_passed]
            input_lon             = self.input_lon[ : idx_time_passed]
            input_spd             = self.input_spd[ : idx_time_passed]
            input_hdg             = self.input_hdg[ : idx_time_passed]
            
            input_time_next       = self.input_time_next[ : idx_time_passed]
            input_lat_next        = self.input_lat_next[ : idx_time_passed]
            input_lon_next        = self.input_lon_next[ : idx_time_passed]
                
            self.check_bird(input_id1, input_bird_size, input_no_inds, input_flock_flag, input_alt)
            #print "id after check", self.id
            # there are probably birds in this set which reached their last timestep.
            # they have to be removed
            # trigger: id1 != id2
            # input_id and self.id have different order - we need the positions
            #  of self.id!
            id_to_remove = input_id1[np.where(input_id1 != input_id2)[0]]
            
            if len(id_to_remove) > 0:
                to_remove = []

                for identity in id_to_remove:
                    
                    
                    if identity in self.id:
                        
                        to_remove = to_remove + list(np.where(self.id == identity)[0])

            # and remove the bird from the simulation
            # but only if there is anything to remove
               # print ("removing from next id")
                self.remove_bird(to_remove) 
                #print "id after remove", self.id

                # input_id1 has to be adjusted as well, otherwise we run into 
                # trouble with the curr_idx-comparison
                # input_id1 = np.delete(input_id1, index_to_remove)


        # position update
        # explanation np.ndenumerate: index is the iterator through the array
        # while id1 is the actual value. 
        # e.g. ([4,3,7]): index = 0,1,2, id1 = 4,3,7
            
            for index, id1 in np.ndenumerate(input_id1):
                index_to_replace = np.where(self.id == id1)[0]        

                self.last_ts[index_to_replace]  = input_time[index]
                self.last_lat[index_to_replace] = input_lat[index]
                self.last_lon[index_to_replace] = input_lon[index]
                
                self.next_ts[index_to_replace]  = input_time_next[index]
                self.next_lat[index_to_replace] = input_lat_next[index]
                self.next_lon[index_to_replace] = input_lon_next[index]                
                
                
                self.tas[index_to_replace]      = input_spd[index]
                self.trk[index_to_replace]      = input_hdg[index]
               # print ('in update, track is', self.trk, 'input hdg', input_hdg)
                        
               # perform some data logging here for the init birdies
                    
                # prepare the data
                #print "birdies are", self.id
                #print "we look at  birdie", input_id1[index], self.id[index_to_replace], type(self.id[index_to_replace]), self.last_lon[index_to_replace], input_lon[index]
                '''bird_data = str(self.id[index_to_replace])  + ' \t ' + str(self.last_lat[index_to_replace]) + \
                    ' \t ' + str(self.last_lon[index_to_replace]) + ' \t ' + str(self.alt[index_to_replace]) + \
                    ' \t ' + str(self.collision_radius[index_to_replace]) + ' \t ' + str(self.bird_size[index_to_replace]) 

                # required: filename,time, bird data: 'date \t time \t id_bird  \t lat \t lon \t alt \t coll_rad \t size \n'
                self.log.write_birds(self.filename2save, simtime, bird_data)
                    
            #save bird trajectories every self.logdt seconds
                if abs(simtime - self.logt0) >= self.logdt:
        
                    self.logt0 = simtime                    
                    
                    # storing only in intervals
                    self.log.save_birds(self.filename2save)'''
                    
            



            self.remove_input(idx_time_passed)
        
        
        self.update_position()
        
        
        return



############################
#CHECK BIRDS
############################

    # check whether the current bird is known already. If not, create    
    def check_bird(self, input_id1, input_bird_size, input_no_inds, input_flock_flag, input_alt):

        # test 1: not in removed_id: if an avian radar bird was eaten by
        # an aircraft, it still might have track data. But as it has been eaten,
        # it can't fly anymore

        # WARNING: THis bird_idx_to add does refer to the idx in the list "input_id1", not the position
        # in the input list containing all birdies
        bird_idx_to_add = np.where(np.in1d(input_id1, self.removed_id, invert = True) & (np.in1d(input_id1, self.id, invert = True)))[0]


#print test[np.where(np.in1d(test, already_in, invert = True) & (np.in1d(test, killed, invert = True)))[0]]

            
            # np.where idx is not in removed or id
            # append the values with these idxs
            
       # print "inputs", input_id1, input_bird_size, input_flock_flag, input_alt
        #print "in removed", np.in1d(input_id1, self.removed_id, invert = True), "in id", (np.in1d(input_id1, self.id, invert = True)), "both", bird_idx_to_add, "id", self.id, "removed", self.removed_id, "input", input_id1      
        
        # collision radius of birds: f(span, size, number)
        add_no_inds = input_no_inds[bird_idx_to_add]
        add_bird_size = input_bird_size[bird_idx_to_add]


        # spans are
        # small: 0.34 m
        # medium: 0.69 m
        # large: 1.43 m
     
        # radius for protected zone around birds
        add_radius = np.zeros(len(bird_idx_to_add))

        # *0.5 because span is diameter and we need radius
        add_radius[add_bird_size == 6] = (np.sqrt(add_no_inds[add_bird_size == 6])* 0.5 * 0.32) + 0.06
        add_radius[add_bird_size == 5] = (np.sqrt(add_no_inds[add_bird_size == 5])* 0.5 * 0.68) + 0.16
        add_radius[add_bird_size == 4] = (np.sqrt(add_no_inds[add_bird_size == 4])* 0.5 * 1.40) + 0.41
        
        add_radius[np.where((add_bird_size == 6)&(add_no_inds ==1))] = 0.5 * 0.32
        add_radius[np.where((add_bird_size == 5)&(add_no_inds ==1))] = 0.5 * 0.68
        add_radius[np.where((add_bird_size == 4)&(add_no_inds ==1))] = 0.5 * 1.40


        
        self.id         = np.append(self.id, input_id1[bird_idx_to_add])
        self.bird_size  = np.append(self.bird_size, add_bird_size )
        self.no_inds    = np.append(self.no_inds, add_no_inds )
        self.flock_flag = np.append(self.flock_flag, input_flock_flag[bird_idx_to_add])
        # and a placeholder for all the other items
        self.last_ts  = np.append(self.last_ts, np.zeros([len(bird_idx_to_add)])) 
        self.last_lat = np.append(self.last_lat, np.zeros([len(bird_idx_to_add)]))
        self.last_lon = np.append(self.last_lon, np.zeros([len(bird_idx_to_add)]))
        
        self.next_ts  = np.append(self.next_ts, np.zeros([len(bird_idx_to_add)])) 
        self.next_lat = np.append(self.next_lat, np.zeros([len(bird_idx_to_add)]))
        self.next_lon = np.append(self.next_lon, np.zeros([len(bird_idx_to_add)]))        
        
        
        
        
        self.lat      = np.append(self.lat, np.zeros([len(bird_idx_to_add)]))
        self.lon      = np.append(self.lon, np.zeros([len(bird_idx_to_add)]))
        self.tas      = np.append(self.tas, np.zeros([len(bird_idx_to_add)]))
        self.trk      = np.append(self.trk, np.zeros([len(bird_idx_to_add)]))

        self.alt      = np.append(self.alt, input_alt[bird_idx_to_add])
        self.collision_radius = np.append(self.collision_radius, add_radius)
        

        return

############
# REMOVE INPUT
############


    def remove_input(self, no_to_remove):

        # remove the info we already looked at
    # these are the first x elements. So the array now starts at the position [element]+1

        self.input_time             = self.input_time[no_to_remove :]
        self.input_id1              = self.input_id1[no_to_remove :]
        self.input_id2              = self.input_id2[no_to_remove :]
        self.input_lat              = self.input_lat[no_to_remove :]
        self.input_lon              = self.input_lon[no_to_remove :]
        self.input_spd              = self.input_spd[no_to_remove :]
        self.input_hdg              = self.input_hdg[no_to_remove :]   
        self.input_alt              = self.input_alt[no_to_remove :]
        self.input_flock_flag       = self.input_flock_flag[no_to_remove :]
        self.input_bird_size        = self.input_bird_size[no_to_remove :]
        self.input_no_inds          = self.input_no_inds[no_to_remove : ]


        self.input_time_next        = self.input_time_next[no_to_remove :]
        self.input_lat_next         = self.input_lat_next[no_to_remove :]
        self.input_lon_next         = self.input_lon_next[no_to_remove :]

        
        return


##############
# UPDATE POSITION
##############

    def update_position(self):
        
        
      
        # timedelta between last and next
        entire_delta_t = self.next_ts - self.last_ts
        delta_time_now = bs.sim.simt - self.last_ts

        self.lat = (self.next_lat - self.last_lat) * (delta_time_now/entire_delta_t) + self.last_lat
        self.lon = (self.next_lon - self.last_lon) * (delta_time_now/entire_delta_t) + self.last_lon
       # print ('last lat', self.last_lat, 'last lon', self.last_lon, 'current lat', self.lat, 'current lon', self.lon, 'next lat', self.next_lat, 'next lon', self.next_lon)
        
       # with open ("taraaaa.txt", "a") as tara:
            
        #    for i in xrange(len(self.id)):
         #       print "writiiiiiing", self.id, self.lat, self.lon, self.last_ts, self.next_ts
          #      tara.write(str(self.id[i]) + ' \t ' + str(simtime) + ' \t ' + str(self.lat[i]) + ' \t ' + str(self.lon[i]) + '\n')
        '''delta_t = abs(simtime - self.last_ts)
        delta_s = self.tas * delta_t
        #print "lastts", self.last_ts[0:10],  "simtime", simtime
        #print "DELTAAAAAAAAAAAAAAAAA", delta_t[0:10]
        
        
        
        self.calculate_position(delta_t, delta_s, simtime)'''
        
        # removal of birds that left the area: only once per minute
        if abs(bs.sim.simt - self.deletet0) >= self.deletedt:
            
            # test if any birds already left the area
            left_south = np.where(self.lat < lat_south)[0]
            left_north = np.where(self.lat > lat_north)[0]
            left_west  = np.where(self.lon < lon_west)[0]
            left_east  = np.where(self.lon > lon_east)[0]
            
            all_left = list(np.unique(np.concatenate((left_south, left_north, left_west, left_east))))
            if len(all_left)  > 0:
              #  print ("removing all left")
                self.remove_bird(all_left)
            
            # reset timer
            self.deletet0 = bs.sim.simt
            
            
        return



######################################################################################################
######################################################################################################      


# CREATE INDIVIDUAL BIRDIES FOR TESTING PURPOSES


######################################################################################################
######################################################################################################        
        
        
        
        
    
    def create_individual(self, birdid, flock_flag, bird_size, no_inds, birdlat, birdlon, birdtrk, birdalt, birdspd):
        '''Creating solo birdies for testing'''
        
        # add one bird object
        n = 1

        # increase number of birds
        self.nbird += n

        # get position of bird
       # print ("bird lat lon input", birdlat, birdlon)
        birdlat = np.array(n * [birdlat])
        birdlon = np.array(n * [birdlon])

        # Limit longitude to [-180.0, 180.0]
        birdlon[birdlon > 180.0] -= 360.0
        birdlon[birdlon < -180.0] += 360.0
        
       # print ('bird lat lon output', birdlat, birdlon)

        # add to birdinfo to lists
        self.id = np.append(self.id, birdid)
        self.flock_flag = np.append(self.flock_flag, flock_flag)
        self.bird_size = np.append(self.bird_size, bird_size)
        self.no_inds = np.append(self.no_inds, no_inds)

        # Positions
        self.lat = np.append(self.lat, birdlat)
        self.lon = np.append(self.lon, birdlon)
        self.alt = np.append(self.alt, birdalt)

        # Heading
        self.trk = np.append(self.trk, birdtrk)

        # Velocities
        self.tas = np.append(self.tas, birdspd)
        
        # TODO: think about vs once we talk 3D birdies
        #vs = 0
        #self.vs = np.append(self.vs, vs)
        
        # radius of protected zone
        if no_inds < 2:

            if bird_size == 6:
                self.collision_radius = np.append(self.collision_radius,  0.5 * 0.32)
            elif bird_size == 5:
                self.collision_radius = np.append(self.collision_radius,  0.5 * 0.68)
            elif bird_size == 4:
                self.collision_radius = np.append(self.collision_radius,  0.5 * 1.40)                
            else:
                print ("solo birdie size not captured", bird_size)
                
                
        else:

            
            if bird_size == 6:
                self.collision_radius = np.append(self.collision_radius,  np.sqrt(no_inds)* 0.5 * 0.32 + 0.06)
            elif bird_size == 5:
                self.collision_radius = np.append(self.collision_radius,  np.sqrt(no_inds)* 0.5 * 0.68 + 0.16)
            elif bird_size == 4:
                self.collision_radius = np.append(self.collision_radius,  np.sqrt(no_inds)* 0.5 * 1.40 + 0.41)                
            else:
                print ("gang birdie size not captured", bird_size)            
            
     
        # and stuff for calculating
        # --> will NOT work yet but since it is for testing only anyways...
        # probably: last lat = input lat, next lat = 10 minute flight with selected heading
        self.last_ts          = np.append(self.last_ts, -999)
        self.last_lat         = np.append(self.last_lat, -999)
        self.last_lon         = np.append(self.last_lon, -999)
        
        self.next_ts          = np.append(self.next_ts, -999)
        self.next_lat         = np.append(self.next_lat, -999)
        self.next_lon         = np.append(self.next_lon, -999)          
                
        
        return

    
    def update_individual_position(self):
        ''' here, the lat lon position will be updated
        
        use the latitude and longitude the bird has right now as well as the current speed and heading.
        
        check the calculate_leo() function on how I did this for scenario birdies. 
        maybe you can even just directly use it - your task to find out :)
            Even if you do, I recommend that you copy the "calculate leo" function, give it another name,
            e.g "calculate_leo_ind_birds" or whatever you like and work with that function. 
            Then you can make adjustments if needed. For example, you can directly access
            self.lat, self.lon, self.trk and do not have to hand them over to the function.
        
        Be careful with radians and degrees and with units!!!
        --> the simulator usually needs degrees and will show imperial units in the user interface
        --> calculations need to be performed in radians and with metric values (see the converting in calculate_leo for example)
        
        
        
        
        '''
        #print ('updating individual position')
        return
    
    
 
    
 
    
    #@core.timed_function(name='example', dt=5)
    #def update(self):
     #       ''' Periodic update function for our example entity. '''
    #        print ('updating, huiiiiiiiii')
            #bird_traf.release_birds()
    
 
    
 
    
    
    
######################################################################################################
######################################################################################################      


# DELETE BIRDS FROM SIMULATION
    # a) because they left the area
    # b) because an aircraft had them for breakfast


######################################################################################################
######################################################################################################      
    
 
    def remove_bird(self, index_to_remove):

        
        # as soon as a bird leaves the simulation, its information has to be removed
        # idx is the index, where the bird info is stored per list
        
        # also gets called when a bird gets hit by an aircraft        

        
        
        
        # list of removed birds
        self.removed_id = np.append(self.removed_id, self.id[index_to_remove])
        
        
        self.nbird = self.nbird - 1 # number of birds
        
        # basic info
        self.id               = np.delete(self.id, index_to_remove)
        self.flock_flag       = np.delete(self.flock_flag, index_to_remove)
        self.bird_size        = np.delete(self.bird_size, index_to_remove)
        self.no_inds          = np.delete(self.no_inds, index_to_remove)
        self.collision_radius = np.delete(self.collision_radius, index_to_remove)

        # Positions
        self.lat              = np.delete(self.lat, index_to_remove)
        self.lon              = np.delete(self.lon, index_to_remove)  
        self.alt              = np.delete(self.alt, index_to_remove)       
        self.trk              = np.delete(self.trk, index_to_remove)  
        

        # Velocities
        self.tas     = np.delete(self.tas, index_to_remove)    # horizontal airspeed [m/s]
       # self.vs     = np.delete(self.vs, index_to_remove)   # vertical speed [m/s]
        
        
        self.last_ts          = np.delete(self.last_ts, index_to_remove)
        self.last_lat         = np.delete(self.last_lat, index_to_remove)
        self.last_lon         = np.delete(self.last_lon, index_to_remove)
        
        self.next_ts          = np.delete(self.next_ts, index_to_remove)
        self.next_lat         = np.delete(self.next_lat, index_to_remove)
        self.next_lon         = np.delete(self.next_lon, index_to_remove)        
        
        return



    
######################################################################################################
######################################################################################################      


# RESET ALL BIRD INFORMATION


######################################################################################################
######################################################################################################  
    
    
    def reset(self):
        # clear all TODO: copy traffarrays
       # print ('BT - in reset')
        self.nbird = 0 # number of birds
        self.scenario_loaded = False
        
        # initialize bid array
        self.id         = np.array([], dtype = int)      
        self.flock_flag = np.array([], dtype = int)
        self.bird_size  = np.array([], dtype = int)
        self.no_inds    = np.array([], dtype = int)
        self.collision_radius = np.array([])
        

        # Positions
        self.lat     = np.array([], dtype=float)  # latitude [deg]
        self.lon     = np.array([], dtype=float)  # longitude [deg]
        self.alt     = np.array([], dtype=float)  # altitude [m]
        self.trk     = np.array([], dtype=float)  # traffic track [deg]

        # Velocities
        self.tas     = np.array([], dtype=float)   # horizontal airspeed [m/s]
        #self.vs     = np.array([], dtype=float)  # vertical speed [m/s]
        
        
         # from old BlueSky
                 # values from file
        self.input_time = np.array([])
        
        # values for calculation 
        self.last_ts    = np.array([])
        self.last_lat   = np.array([])
        self.last_lon   = np.array([])
        
        self.next_ts    = np.array([])
        self.next_lat   = np.array([])
        self.next_lon   = np.array([])        
        
        self.input_time_next = np.array([])
        self.input_lat_next = np.array([])
        self.input_lon_next = np.array([])
        


        
        # birds which experienced a strike don't fly anymore
        self.removed_id = np.array([])
        
        # set timer for deletions
        self.deletedt = 60. # [s] interval to delete birds that left the area
        self.deletet0 = -self.deletedt # init


        # scheduling of datalog
        self.logdt = 10.
        self.logt0 = -self.logdt
        
        
     #   print ("BT - END RESET", self.id)
        
        
        return
    
######################################################################################################
######################################################################################################      


# HELPER FUNCTIONS


######################################################################################################
######################################################################################################   
    



######################################################################################################      
# id2idx
######################################################################################################
    
    def id2idx(self, birdid):
        """Find index of bird id"""

        return np.where(self.id == np.char.upper(np.array(birdid)))[0][0]




######################################################################################################      
# assigning input values for bird scenarios
######################################################################################################
    def assign_values(self, data):
        self.input_id1 = np.array(pd.to_numeric(data["id"])).astype(int)
        self.input_id2 = np.array(pd.to_numeric(data["id1"])).astype(int)
        self.input_lat = np.array(pd.to_numeric(data["lat"]))
        self.input_lon = np.array(pd.to_numeric(data["lon"]))
        self.input_spd = np.array(pd.to_numeric(data["spd"]))
        self.input_hdg = np.array(pd.to_numeric(data["hdg"]))
        self.input_alt = np.array(pd.to_numeric(data["alt"]))
        
        self.input_bird_size  = np.array(pd.to_numeric(data["cat"])).astype(int)
        self.input_flock_flag = np.array(data["flock_flag"])
        self.input_no_inds    = np.array(data['no_individuals']).astype(int)
        
        self.input_time = np.array(pd.to_numeric(data['date']))   
        
        
        self.input_time_next = np.array(pd.to_numeric(data['timeshift']))
        self.input_lat_next = np.array(pd.to_numeric(data['latshift']))
        self.input_lon_next = np.array(pd.to_numeric(data['lonshift']))
        
        

        return

######################################################################################################      
# Leo's positions
###################################################################################################### 


    def calculate_leo_position(self, delta_t, delta_s, lat_in, lon_in, hdg):
        '''to prevent birdies from being deleted one timestep too early respective start flying in a weird direction,
        we give them a theoretical final position after which they are deleted '''
        theta = np.radians(hdg)
        last_lat = np.radians(lat_in)
        last_lon = np.radians(lon_in)
       # print ('in leo', lat_in, lon_in, hdg)
        
        d = delta_s/earth_radius


        # calculate the next position with the haversine function
        # check http://www.movable-type.co.uk/scripts/latlong.html for reference
        
        lat_pos = np.arcsin(np.sin(last_lat)*np.cos(d) + np.cos(last_lat)*np.sin(d)*np.cos(theta))
        
        lon_pos = last_lon + np.arctan2(np.sin(theta)*np.sin(d)*np.cos(last_lat),
                      np.cos(d) - np.sin(last_lat)*np.sin(lat_pos))        
                
        lat_theo = np.degrees(lat_pos)
        lon_theo = np.degrees(lon_pos)
        #print lat_theo
        #with open ("taraaaa.txt", "a") as tara:
            

        return lat_theo, lon_theo



    
    
######################################################################################################      
# sending birds from sim to gui
######################################################################################################    
    def release_birds(self):
         
         '''release them to the visual world '''
        
         ''' we need to differentiate between individuals and scenarios
        - for now, either one or the other'''
         
         if self.scenario_loaded:
            
            self.update_scen()
            
         else:
            '''by now, individual birdies do not move
            here is the call to do updates of individual bird positions'''
            self.update_individual_position()
            
        
         data = dict()
         # id is necessary for some gui stuff
         data['id']         = self.id
         #data['type']       = self.type
         data['lat']        = self.lat
         data['lon']        = self.lon
         data['alt']        = self.alt
         data['trk']        = self.trk
         #data['tas']         = self.tas
         #data['hs']         = self.hs
       #  print ('bird track is', self.trk)
         
         
 
         # send bird data
         bs.net.send_stream(b'BIRDDATA', data)
        # print ("BT - releasing", data['id'])
        
         return   



