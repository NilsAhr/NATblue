"""
plugin with update function in 3h intervals
ERA5 data was migrated to CDS completely https://cds.climate.copernicus.eu

- you need an CDS / ECMWF account to retrieve data!!!
- usage: [windecmwfup 20 0 90 -90 yyyy mm dd hh]
        NAT [10 50 90 -130 yyyy mm dd hh]
- resolution of data need to be 0.25x0.25 (see reshapefactor)

written by Nils Ahrenhold (TUD/DLR) 22.01.2025 """

from pathlib import Path
import cdsapi
import datetime
import numpy as np
import bluesky as bs
import netCDF4 as nc
from bluesky import stack
from bluesky.core import timed_function
from bluesky.traffic.windsim import WindSim


datadir = Path('')


def init_plugin():
    global datadir
    datadir = bs.resource(bs.settings.data_path) / 'NetCDF'

    if not datadir.is_dir():
        datadir.mkdir()

    global windecmwfup
    windecmwfup = WindECMWFUP()

    config = {
        'plugin_name': 'WINDECMWFUP',
        'plugin_type': 'sim'
    }

    return config

class WindECMWFUP(WindSim):
    def __init__(self):
        super().__init__()
        self.year  = 0
        self.month = 0
        self.day   = 0
        self.hour  = 0
        self.lat0  = -90 # South Pole (min latitude)
        self.lon0  = -180 # Western Hemisphere (min longitude)
        self.lat1  = 90 # North Pole (max latitude)
        self.lon1  = 180 # Eastern Hemisphere (max longitude)

        # Switch for periodic loading of new GFS data
        self.autoload = False
        
    def fetch_nc(self, year, month, day):
        """
        Retrieve weather data via the CDS API for multiple pressure levels
        """
        
        ymd = "%04d%02d%02d" % (year, month, day)
        fname = f'p_levels_{ymd}.nc'
        fpath = datadir / fname

        # Use cached file if already loaded for this date
        if hasattr(self, "netcdf_date") and self.netcdf_date == ymd:
            return self.netcdf  # Use already loaded file
    
        # Close old file if open and switching to a new date
        if hasattr(self, "netcdf") and hasattr(self, "netcdf_date") and self.netcdf is not None:
            try:
                self.netcdf.close()
            except Exception:
                pass    
        # **Cache file in memory** instead of reading from disk each time
        #if hasattr(self, "netcdf") and self.netcdf.getncattr("date") == ymd:
        #    return self.netcdf  # Use already loaded file

        print(f"Looking for NetCDF file at: {fpath}")
        if not fpath.is_file():
            bs.scr.echo("Downloading file, please wait...")
    
            # Set client
            c = cdsapi.Client()
            
            # Retrieve data 
            c.retrieve(
                'reanalysis-era5-pressure-levels',
                {
                    'product_type': ['reanalysis'],
                    'variable': [
                        'u_component_of_wind',
                        'v_component_of_wind'
                    ],
                    'year': year,
                    'month': month,
                    'day': day,
                    'time': [
                        '00:00', '01:00', '02:00', '03:00',
                        '04:00', '05:00', '06:00', '07:00',
                        '08:00', '09:00', '10:00', '11:00',
                        '12:00', '13:00', '14:00', '15:00',
                        '16:00', '17:00', '18:00', '19:00',
                        '20:00', '21:00', '22:00', '23:00'
                    ],
                    'pressure_level': [
                        '100', '125', '150', 
                        '175', '200', '225',
                        '250', '300', '350',
                        '400', '450', '500',
                        '550', '600', '650', 
                        '700', '750', '775',
                        '800'
                    ],
                    'data_format': 'netcdf',
                    "download_format": "unarchived",
                    "area": [90, -180, -90, 180]             
                },
                fpath)
    
        bs.scr.echo("Download completed.")
        netcdf = nc.Dataset(fpath, mode='r')
        self.netcdf = netcdf
        self.netcdf_date = ymd  # <-- Store the date string
        return netcdf

    
    def extract_wind(self, netcdf, lat0, lon0, lat1, lon1, hour):

        # Load reanalysis data 
        level = netcdf['pressure_level'][:].data
        lats  = netcdf['latitude'][:].data
        lons  = netcdf['longitude'][:].data
        vxs_  = netcdf['u'][:].squeeze().data
        vys_  = netcdf['v'][:].squeeze().data
        
        # Close data for performance
        #netcdf.close()   
        
        # Transform pressure levels to altitude
        p = level * 100
        h = (1 - (p / 101325.0)**0.190264) * 44330.76923    # in meters
        
        # Set hour to rounded hour
        #hour = round(hour/3)
        
        # Construct 2D array of all data points
        lats_ = np.tile(np.repeat(lats, len(lons)), len(level))
        lons_ = np.tile(lons, len(lats)*len(level))
        alts_ = np.repeat(h, len(lats)*len(lons))       
        vxs_  = vxs_[hour,:,:,:].flatten() #u_component_of_wind - selects wind data for a specific time step (given by hour)
        vys_  = vys_[hour,:,:,:].flatten() #v_component_of_wind - selects wind data for a specific time step (given by hour)
            
        # Convert longitudes
        lons_ = (lons_ + 180) % 360.0 - 180.0     # convert range from 0~360 to -180~180

        # Reduce area based on lat lon limits
        lat0_ = min(lat0, lat1)
        lat1_ = max(lat0, lat1)
        lon0_ = min(lon0, lon1)
        lon1_ = max(lon0, lon1)

        mask = (lats_ >= lat0_) & (lats_ <= lat1_) & (lons_ >= lon0_) & (lons_ <= lon1_)

        data = np.array([lats_[mask], lons_[mask], alts_[mask], vxs_[mask], vys_[mask]])

        return data


    @stack.command(name='WINDECMWFUP')
    def loadwind(self, lat0: 'lat', lon0: 'lon', lat1: 'lat', lon1: 'lon',
               year: int=None, month: int=None, day: int=None, hour: int=None):
        ''' WINDECMWFUP: Load a windfield directly from CDS database.


            Arguments:
            - lat0 (south), lon0(east), lat1(north), lon1(west) [deg]:
            - windecmwfup 20 0 90 -90 yyyy mm dd hh

            Bounding box in which to generate wind field
            - year, month, day, hour: Date and time of wind data (optional, will use
              current simulation UTC if not specified).
        '''
        self.lat0, self.lon0, self.lat1, self.lon1 =  min(lat0, lat1), \
                              min(lon0, lon1), max(lat0, lat1), max(lon0, lon1)
        self.year = year or bs.sim.utc.year
        self.month = month or bs.sim.utc.month
        self.day = day or bs.sim.utc.day
        #self.hour = hour or bs.sim.utc.hour
        self.hour = hour if hour is not None else bs.sim.utc.hour  # <-- Only override if hour is not provided

        # round hour to 3 hours
        # self.hour  = round(self.hour/3) * 3
        
        if self.hour == 24:
            ymd0 = "%04d%02d%02d" % (self.year, self.month, self.day)
            ymd1 = (datetime.datetime.strptime(ymd0, '%Y%m%d') + 
                    datetime.timedelta(days=1))
            self.year  = ymd1.year
            self.month = ymd1.month
            self.day   = ymd1.day
            self.hour  = 0

        txt = "Loading wind field for %s-%s-%s-%02d:00..." % (self.year, self.month, self.day, self.hour)
        bs.scr.echo("%s" % txt)

        netcdf = self.fetch_nc(self.year, self.month, self.day)

        if netcdf is None or self.lat0 == self.lat1 or self.lon0 == self.lon1:
            return False, "Wind data non-existend in area [%d, %d], [%d, %d]. " \
                % (self.lat0, self.lat1, self.lon0, self.lon1) \
                + "time: %04d-%02d-%02d" \
                % (self.year, self.month, self.day)

        # first clear exisiting wind field
        self.clear()

        # add new wind field
        data = self.extract_wind(netcdf, self.lat0, self.lon0, self.lat1, self.lon1, self.hour).T
        
        data = data[np.lexsort((data[:, 2], data[:, 1], data[:, 0]))] # Sort by lat, lon, alt        
        reshapefactor = int((1 + (max(self.lat0, self.lat1) - min(self.lat0, self.lat1))*4) * \
                            (1 + (max(self.lon0, self.lon1) - min(self.lon0, self.lon1))*4))
        
        lat     = np.reshape(data[:,0], (reshapefactor, -1)).T[0,:]
        lon     = np.reshape(data[:,1], (reshapefactor, -1)).T[0,:]
        veast   = np.reshape(data[:,3], (reshapefactor, -1)).T
        vnorth  = np.reshape(data[:,4], (reshapefactor, -1)).T
        windalt = np.reshape(data[:,2], (reshapefactor, -1)).T[:,0]

        self.addpointvne(lat, lon, vnorth, veast, windalt)        

        self.autoload = True  # Enable autoload for next update
        return True, "Wind field updated in area [%d, %d], [%d, %d]. " \
            % (self.lat0, self.lat1, self.lon0, self.lon1) \
            + "time: %04d-%02d-%02d-%02d:00" \
            % (self.year, self.month, self.day, self.hour)

    @timed_function(name='WINDECMWFUP', dt=3600) #1h = 3600, 2h = 7200, 3h = 10800,  4h = 14400, 5h = 18000, 6h = 21600, 8h = 28800
    def update(self):
        if self.autoload:
            bs.scr.echo("updating windfield")
            # Increment the hour by to get the next timestep
            self.hour += 1
            
            # Check if the hour exceeds 24, and adjust the date if needed
            if self.hour >= 24:
                self.hour = 0
                current_date = datetime.date(self.year, self.month, self.day)
                next_date = current_date + datetime.timedelta(days=1)
                self.year = next_date.year
                self.month = next_date.month
                self.day = next_date.day

            # **Only reload data if date changed**
            if self.hour == 0:
                self.netcdf = self.fetch_nc(self.year, self.month, self.day)

            info = f"Current wind data hour: {self.hour:02d}:00"
            bs.scr.echo(info)
            # Load the wind data for the next time step for defined coordinates
            _, txt = self.loadwind(self.lat0, self.lon0, self.lat1, self.lon1,
                                   self.year, self.month, self.day, self.hour)

            bs.scr.echo("%s" % txt)
            bs.scr.echo(f"Wind field updated: {self.year}-{self.month}-{self.day} {self.hour}:00")
            
            """
            # **Reuse extracted wind data**
            wind_data = self.extract_wind(self.netcdf, self.lat0, self.lon0, self.lat1, self.lon1, self.hour)

            # First clear existing wind field
            self.clear()

            # Reshape and apply data
            lat, lon, alt, veast, vnorth = wind_data.T
            self.addpointvne(lat, lon, vnorth, veast, alt)
            """
                        
            
        