# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 15:05:13 2017

@author: Isabel
"""
import numpy as np
import pandas as pd

def randomize_birds(seed, filename):

    
    # first get the file
    #df = pd.read_csv('bird_movement_plan.csv', sep = "\t", \
    #                     names = ['id', 'time', 'lon','lat', 'alt', 'size', 'number', 'ff',\
    #                     'id1', 'spd', 'hdg'], index_col = False)
    
    df = pd.read_csv(filename, sep = "\,", \
                         names = ['id', 'date', 'lon','lat', 'alt', 'cat', 'no_individuals', 'flock_flag',\
                         'id1', 'hdg', 'spd', 'lat_s1', 'lat_n1', 'lon_w1', 'lon_e1','lat_s2', 'lat_n2', 'lon_w2', 'lon_e2'], index_col = False, engine = 'python')


    #print "at init", df['id']
    # place the seed
    np.random.seed(seed)
    
    # make sure that we only use the first occurrence of every birdie :)
    df['duplicate'] =  df.duplicated(['id'])
    
    # df1 id the one with the first occurrences
    df1 = df[df.duplicated(['id']) == False]
    # df2 are the ones to follow
    df2 = df[df.duplicated(['id']) == True]

    
    # altitude
        
        # for all birdies
        
        # BUT ONLY FOR FIRST OCCURRENCE OF EVERY BIRD
    
        
        # respective altitude band --> check if the inclusive is at lower or upper boundary
        
        # 0-200
    mask = (df1['alt'] > 0.) & (df1['alt'] <= 200.)
    
    
    
    df1.loc[mask, 'alt'] = np.random.uniform(0., 200., len(df1.loc[mask, 'alt']))
    
    
        
        # 200-400
    mask = (df1['alt'] > 200.) & (df1['alt'] <= 400.)
    
    df1.loc[mask, 'alt'] = np.random.uniform(200., 400., len(df1.loc[mask, 'alt']))    
        # 400-600
    mask = (df1['alt'] > 400.) & (df1['alt'] <= 600.)
    
    df1.loc[mask, 'alt'] = np.random.uniform(400., 600., len(df1.loc[mask, 'alt']))       
        # 600-800
    mask = (df1['alt'] > 600.) & (df1['alt'] <= 800.)
    
    df1.loc[mask, 'alt'] = np.random.uniform(600., 800., len(df1.loc[mask, 'alt']))       
        # 800-1000
    mask = (df1['alt'] > 800.) & (df1['alt'] <= 1000.)
    
    df1.loc[mask, 'alt'] = np.random.uniform(800., 1000., len(df1.loc[mask, 'alt']))   
    # latitude and longitude
    
        #   only for weather radar birdies (alt >=200m)
    
        #    creates in the field: lat_south, lat_north, lon_west, lon_east
    
        #   fly-in birdies: within their field... Add parameters? Then we need to add them everywhere (or is empty)
    
    
    # speed: only weather radar birdies (alt >=200m) 
    
            #  plusminus 6 around input value
    
    df1['spd_low'] = df1['spd'] - 6.
    # careful: spd_low can become negative. Where this happens, spd_low is replaced with spd
    mask = df1['spd_low'] <= 0.

    df1.loc[mask, 'spd_low'] = df1['spd']
    df1['spd_high'] = df1['spd'] + 6.
    

    #df1.loc['alt' > 200., 'spd'] = np.random.uniform(df1['spd_low'], df1['spd_high'], len(df1.loc['alt' > 200., 'spd']))
    #df1.loc[(df1['alt'] > 200.),'spd_new'] = np.random.uniform(df1['spd_low'], df1['spd_high'], len(df1.loc[df1['alt'] > 200.]))
    
    # HAS TO BE spd EVENTUALLY!!!!
    df1['spd_new'] = np.random.uniform(df1['spd_low'], df1['spd_high'], len(df1))
    df1.loc[(df1['alt'] > 200.),'spd'] = df1['spd_new']
    
    # remove non-needed columns
    del df1['spd_new']
    del df1['spd_low']
    del df1['spd_high']
    
    
    # direction: only weather radar birdies (alt >=200m) 
    
            # TAKE CARE OF 360 degrees!!!
            # plusminus 45 degrees from init value
            # standard deviation for heading is plusminus 45 degrees
    df1.loc[(df1['alt'] > 200.), 'random_hdg'] = np.random.uniform(-45., 45., len(df1.loc[df1['alt'] > 200.]))       
    df1.loc[(df1['alt'] > 200.), 'hdg'] = (df1['random_hdg'] + df1['hdg'] + 360.) % 360
    
    
    
    df1 = df1.drop(['random_hdg'], axis = 1)
    
    
    # latitude and longitude: only for weather radar birdies
    #  first: set minima and maxima for lat and lon. 
    # Random choice because for birds flying in diagonally (e.g from North-West)
    # can do so from two directions (North or West in this case)
    # it is a 50:50 chance from which side they come.
    # birds from one side or within the area: 1- and 2-values are identical
    
    # pairwise needed... either it is lat_s1, lat_n1, lon_w1, lon_e1 OR lat_s2, lat_n2, lon_w2, lon_e2
    df1['designator'] = np.random.choice([1,2], len(df1))
    
    # pair 1
    df1.loc[(df1['designator'] == 1), 'lat'] =np.random.uniform(df1.loc[(df1['designator'] == 1), 'lat_s1'], df1.loc[(df1['designator'] == 1), 'lat_n1'], len(df1.loc[df1['designator'] == 1]))
    df1.loc[(df1['designator'] == 1), 'lon'] = np.random.uniform(df1.loc[(df1['designator'] == 1), 'lon_w1'], df1.loc[(df1['designator'] == 1), 'lon_e1'], len(df1.loc[df1['designator'] == 1]))
    
    # pair 2
    df1.loc[(df1['designator'] == 2), 'lat'] =np.random.uniform(df1.loc[(df1['designator'] == 2), 'lat_s2'], df1.loc[(df1['designator'] == 2), 'lat_n2'], len(df1.loc[df1['designator'] == 2]))
    df1.loc[(df1['designator'] == 2), 'lon'] = np.random.uniform(df1.loc[(df1['designator'] == 2), 'lon_w2'], df1.loc[(df1['designator'] == 2), 'lon_e2'], len(df1.loc[df1['designator'] == 2]))
    
    
    
    
    
    
    
    # remove non-needed columns
    df1 = df1.drop(['lat_s1', 'lat_s2', 'lat_n1', 'lat_n2', 'lon_w1', 'lon_w2', 'lon_e1', 'lon_e2','duplicate', 'designator'], axis = 1)
    df2 = df2.drop(['lat_s1', 'lat_s2', 'lat_n1', 'lat_n2', 'lon_w1', 'lon_w2', 'lon_e1', 'lon_e2', 'duplicate'], axis = 1)
    
  #  print "before merging",df1['id']
    
    
    
    # merge the frames back together
    
    df_merged = pd.concat([df1, df2])

    # for the "next_ts", lat, lon, alt, we need to shift the respective cells
    # BUT FIRST SORT CORRECTLY!
    df_merged = df_merged.sort_values(['id','date' ], ascending=[True, True])

    df_merged['timeshift'] = df_merged['date'].shift(-1)
    df_merged['latshift'] = df_merged['lat'].shift(-1)
    df_merged['lonshift'] = df_merged['lon'].shift(-1)    


    #why not directly convert?

    df_merged['lat']      = np.radians(df_merged['lat'])
    df_merged['lon']      = np.radians(df_merged['lon'])
    df_merged['latshift'] = np.radians(df_merged['latshift'])
    df_merged['lonshift'] = np.radians(df_merged['lonshift'])    
    
    df_merged = df_merged.sort_values(by='date') 
   # print "end of randomizing ", df_merged['id']
    
    return df_merged
    
    
#df = randomize_birds(42, "bla")
#print df
