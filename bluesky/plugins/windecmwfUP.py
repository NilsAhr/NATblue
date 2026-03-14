"""
plugin with update function in configurable intervals
ERA5 data was migrated to CDS completely https://cds.climate.copernicus.eu

- you need an CDS / ECMWF account to retrieve data!!!
- usage: [windecmwfup 20 0 90 -90 yyyy mm dd hh]
        NAT [10 50 90 -130 yyyy mm dd hh]
        with resolution: [windecmwfup 20 0 90 -90 yyyy mm dd hh 1.0 10800]
- native resolution of data is 0.25x0.25 (see reshapefactor)
- optional spatial_res (deg) and time_res (s) control subsampling
  and update rate; defaults are 0.25 deg and 3600 s (1h)

written by Nils Ahrenhold (TUD/DLR) 22.01.2025 """

from pathlib import Path
import cdsapi
import datetime
import numpy as np
import bluesky as bs
import netCDF4 as nc
from bluesky import stack
from bluesky.core import timed_function
from bluesky.core.simtime import TimerMeta
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

        # Resolution settings (defaults = native 0.25 deg, 1h update)
        self.spatial_res = 0.25   # spatial resolution in degrees
        self.time_res    = 3600   # update interval in seconds

        # Switch for periodic loading of new GFS data
        self.autoload = False
        #self.just_loaded = False
        
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
            stack.echo("Downloading file, please wait...")
    
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
                        '1000','975','950','925','900','875','850','825','800','775',
                        '750','700','650','600','550','500','450','400','350','300',
                        '275','250','225','200','175','150','125','100'
                    ],
                    # 'pressure_level': [
                    #     '100', '125', '150', 
                    #     '175', '200', '225',
                    #     '250', '300', '350',
                    #     '400', '450', '500',
                    #     '550', '600', '650', 
                    #     '700', '750', '775',
                    #     '800'
                    # ],
                    'data_format': 'netcdf',
                    "download_format": "unarchived",
                    "area": [90, -180, -90, 180]             
                },
                fpath)
    
        stack.echo("Download completed.")
        netcdf = nc.Dataset(fpath, mode='r')
        self.netcdf = netcdf
        self.netcdf_date = ymd  # <-- Store the date string
        return netcdf

    
    def extract_wind(self, netcdf, lat0, lon0, lat1, lon1, hour, spatial_res=0.25):

        # Load reanalysis data 
        level = netcdf['pressure_level'][:].data
        lats  = netcdf['latitude'][:].data
        lons  = netcdf['longitude'][:].data
        vxs_  = netcdf['u'][:].squeeze().data
        vys_  = netcdf['v'][:].squeeze().data

        # Subsample lat/lon if spatial_res > native 0.25 deg
        native_res = 0.25
        step = max(1, int(round(spatial_res / native_res)))
        if step > 1:
            lats = lats[::step]
            lons = lons[::step]
            vxs_ = vxs_[:, :, ::step, ::step]   # dims: [time, level, lat, lon]
            vys_ = vys_[:, :, ::step, ::step]
        
        # Transform pressure levels to altitude
        p = level * 100
        h = (1 - (p / 101325.0)**0.190264) * 44330.76923    # in meters
        
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
               year: int=None, month: int=None, day: int=None, hour: int=None,
               spatial_res: float=None, time_res: int=None):
        ''' WINDECMWFUP: Load a windfield directly from CDS database.


            Arguments:
            - lat0, lon0, lat1, lon1 [deg]: [two corner points]
            - windecmwfup 20 0 90 -90 yyyy mm dd hh
            - windecmwfup 20 0 90 -90 yyyy mm dd hh spatial_res time_res
              e.g. windecmwfup 10 50 90 -130 2025 03 14 08 1.0 10800

            Bounding box in which to generate wind field
            - year, month, day, hour: Date and time of wind data (optional, will use
              current simulation UTC if not specified).
            - spatial_res [deg]: Spatial resolution for wind field subsampling
              (default 0.25 = native ERA5 resolution). E.g. 1.0 for 1 deg.
            - time_res [seconds]: Wind field update interval in seconds
              (default 3600 = 1h). E.g. 10800 for 3h updates.
        '''
        self.lat0, self.lon0, self.lat1, self.lon1 =  min(lat0, lat1), \
                              min(lon0, lon1), max(lat0, lat1), max(lon0, lon1)
        self.year = year or bs.sim.utc.year
        self.month = month or bs.sim.utc.month
        self.day = day or bs.sim.utc.day
        #self.hour = hour or bs.sim.utc.hour
        self.hour = hour if hour is not None else bs.sim.utc.hour  # <-- Only override if hour is not provided

        # Apply spatial and time resolution if provided
        if spatial_res is not None:
            self.spatial_res = spatial_res
        if time_res is not None:
            self.time_res = int(time_res)
            # Dynamically adjust the timed_function update interval
            timer = TimerMeta.gettimer('WINDECMWFUP')
            if timer is not None:
                timer.setdt(self.time_res)
                stack.echo(f"Wind update interval set to {self.time_res}s "
                           f"({self.time_res/3600:.1f}h)")

        stack.echo(f"Resolution: spatial={self.spatial_res} deg, "
                   f"update interval={self.time_res}s ({self.time_res/3600:.1f}h)")

        if self.hour == 24:
            ymd0 = "%04d%02d%02d" % (self.year, self.month, self.day)
            ymd1 = (datetime.datetime.strptime(ymd0, '%Y%m%d') + 
                    datetime.timedelta(days=1))
            self.year  = ymd1.year
            self.month = ymd1.month
            self.day   = ymd1.day
            self.hour  = 0

        txt = "Loading wind field for %s-%s-%s-%02d:00..." % (self.year, self.month, self.day, self.hour)
        stack.echo("%s" % txt)

        netcdf = self.fetch_nc(self.year, self.month, self.day)

        if netcdf is None or self.lat0 == self.lat1 or self.lon0 == self.lon1:
            return False, "Wind data non-existend in area [%d, %d], [%d, %d]. " \
                % (self.lat0, self.lat1, self.lon0, self.lon1) \
                + "time: %04d-%02d-%02d" \
                % (self.year, self.month, self.day)

        # first clear exisiting wind field
        self.clear()

        # add new wind field
        data = self.extract_wind(netcdf, self.lat0, self.lon0, self.lat1, self.lon1,
                                 self.hour, spatial_res=self.spatial_res).T
        
        data = data[np.lexsort((data[:, 2], data[:, 1], data[:, 0]))] # Sort by lat, lon, alt
        # Compute reshapefactor from actual resolution (points per degree = 1/spatial_res)
        ppd = 1.0 / self.spatial_res  # points per degree
        reshapefactor = int((1 + (max(self.lat0, self.lat1) - min(self.lat0, self.lat1)) * ppd) * \
                            (1 + (max(self.lon0, self.lon1) - min(self.lon0, self.lon1)) * ppd))
        
        lat     = np.reshape(data[:,0], (reshapefactor, -1)).T[0,:]
        lon     = np.reshape(data[:,1], (reshapefactor, -1)).T[0,:]
        veast   = np.reshape(data[:,3], (reshapefactor, -1)).T
        vnorth  = np.reshape(data[:,4], (reshapefactor, -1)).T
        windalt = np.reshape(data[:,2], (reshapefactor, -1)).T[:,0]

        self.addpointvne(lat, lon, vnorth, veast, windalt)        

        #self.just_loaded = True  # Mark that data was just loaded
        self.autoload = True  # Enable autoload for next update
        return True, "Wind field updated in area [%d, %d], [%d, %d]. " \
            % (self.lat0, self.lat1, self.lon0, self.lon1) \
            + "time: %04d-%02d-%02d-%02d:00 (res: %.2f deg, %ds update)" \
            % (self.year, self.month, self.day, self.hour, self.spatial_res, self.time_res)

    @timed_function(name='WINDECMWFUP', dt=3600) #default 1h; changed dynamically via time_res parameter
    def update(self):
        if self.autoload:
            # Advance by the number of hours corresponding to time_res,
            # anchored to the initial load time (not synced to sim UTC,
            # which avoids misalignment with the optimiser's time grid)
            hours_step = max(1, self.time_res // 3600)
            self.hour += hours_step

            # Handle day rollover
            if self.hour >= 24:
                days_ahead = self.hour // 24
                self.hour = self.hour % 24
                current_date = datetime.date(self.year, self.month, self.day)
                next_date = current_date + datetime.timedelta(days=days_ahead)
                self.year = next_date.year
                self.month = next_date.month
                self.day = next_date.day
            
            # Only reload NetCDF if the date changed    
            ymd_now = f"{self.year:04d}{self.month:02d}{self.day:02d}"
            if not hasattr(self, "netcdf_date") or self.netcdf_date != ymd_now:
                self.netcdf = self.fetch_nc(self.year, self.month, self.day)

            stack.echo(f"Updating wind field to {self.year}-{self.month:02d}-{self.day:02d} {self.hour:02d}:00")
            _, txt = self.loadwind(self.lat0, self.lon0, self.lat1, self.lon1,
                                self.year, self.month, self.day, self.hour)
            stack.echo(txt)

            # Check if the hour exceeds 24, and adjust the date if needed
            #if self.hour >= 24:
            #    self.hour = 0
            #    current_date = datetime.date(self.year, self.month, self.day)
            #    next_date = current_date + datetime.timedelta(days=1)
            #    self.year = next_date.year
            #    self.month = next_date.month
            #    self.day = next_date.day

            # **Only reload data if date changed**
            #if self.hour == 0:
            #    self.netcdf = self.fetch_nc(self.year, self.month, self.day)

            #info = f"Current wind data hour: {self.hour:02d}:00"
            #stack.echo(info)
            # Load the wind data for the next time step for defined coordinates
            #_, txt = self.loadwind(self.lat0, self.lon0, self.lat1, self.lon1,
            #                       self.year, self.month, self.day, self.hour)

            #stack.echo("%s" % txt)
            #stack.echo(f"Wind field updated: {self.year}-{self.month}-{self.day} {self.hour}:00")
            
            """
            # **Reuse extracted wind data**
            wind_data = self.extract_wind(self.netcdf, self.lat0, self.lon0, self.lat1, self.lon1, self.hour)

            # First clear existing wind field
            self.clear()

            # Reshape and apply data
            lat, lon, alt, veast, vnorth = wind_data.T
            self.addpointvne(lat, lon, vnorth, veast, alt)
            """
                        
            
        