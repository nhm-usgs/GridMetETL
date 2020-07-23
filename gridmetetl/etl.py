"""Main module."""
import geopandas
import pandas as pd
from netCDF4 import default_fillvals, Dataset
from numpy import arange, dtype, float32, zeros, asarray
import sys
import xarray as xr
from helper import get_gm_url, np_get_wval, getaverage
import requests
from requests.exceptions import HTTPError
from datetime import datetime
from pathlib import Path
import numpy as np
import netCDF4


class FpoNHM:
    """ Class for fetching climate data and parsing into netcdf
        input files for use with the USGS operational National Hydrologic
        Model (oNHM).  Workflow:
            1) Initialize(): fetch climate data
            2) Run(): map/interpolate onto hru
            3) Finalize(): write netcdf input files
        Mapping options:
            1) weighted average based on intersection area of hru
                with netcdf file cells.
            2) rasterstats - zonal averaging

    """

    def __init__(self, climsource='GridMetSS'):
        """
        Initialize class

        :param  : number of days past to retrieve
        :param climsource: Constant for now but may have multiple
            choice for data sources in the future.  Currently default is
            GridMet:  http://www.climatologylab.org/gridmet.html
        """
        self.wghts_id = None
        self.climsource = climsource
        if climsource == 'GridMetSS':
            self.gmss_vars = {
                'tmax': 'daily_maximum_temperature',
                'tmin': 'daily_minimum_temperature',
                'ppt': 'precipitation_amount',
                'rhmax': 'daily_maximum_relative_humidity',
                'rhmin': 'daily_minimum_relative_humidity',
                'ws': 'daily_mean_wind_speed',
                'srad': 'daily_mean_shortwave_radiation_at_surface'}
        self.partial = False
        self.vars = ['tmax', 'tmin', 'ppt', 'rhmax', 'rhmin', 'ws']

        # type of retrieval (days) retrieve by previous number of days - used in operational mode
        # or (date) used to retrieve specific period of time
        self.type = None

        self.numdays = None

        # prefix for file names - default is ''.
        self.fileprefix = None

        # xarray containers for tempurature max, temperature min and precipitation
        self.dstmax = None
        self.dstmin = None
        self.dsppt = None
        self.dsrhmax = None
        self.dsrhmin = None
        self.dsws = None
        self.dssrad = None
        self.ds = None

        # Geopandas dataframe that will hold hru id and geometry
        self.gdf = None
        self.gdf1 = None

        # input and output path directories
        self.iptpath = None
        self.optpath = None

        # weights file
        self.wghts_file = None

        # start and end dates of using type == date
        self.start_date = None
        self.end_date = None

        # handles to netcdf climate data
        # Coordinates
        self.lat_h = None
        self.lon_h = None
        self.time_h = None
        # Geotransform
        self.crs_h = None
        # Climate data
        self.tmax_h = None
        self.tmin_h = None
        self.tppt_h = None
        self.rhmax_h = None
        self.rhmin_h = None
        self.ws_h = None
        # Dimensions
        self.dayshape = None
        self.lonshape = None
        self.latshape = None

        # num HRUs
        self.num_hru = None

        # grouby hru_id_nat on wieghts file
        self.unique_hru_ids = None

        # numpy arrays to store mapped climate data
        self.np_tmax = None
        self.np_tmin = None
        self.np_ppt = None
        self.np_rhmax = None
        self.np_rhmin = None
        self.np_ws = None
        self.np_srad = None

        # logical use_date
        self.use_date = False

        # Starting date based on numdays
        self.str_start = None

        self.dstmax = None
        self.dstmin = None
        self.dsrhmax = None
        self.dsrhmin = None
        self.dsws = None
        self.dssrad = None
        self.dsppt = None

    def write_extract_file(self, ivar, incfile, url, params):
        file = requests.get(url, params=params)
        file.raise_for_status()
        tfile = self.iptpath / (self.fileprefix + ivar + (datetime.now().strftime('%Y_%m_%d')) + '.nc')
        incfile.append(tfile)
        with open(tfile, 'wb') as fh:
            fh.write(file.content)
        fh.close()

    def initialize(self, partial, ivars, iptpath, optpath, weights_file, etype=None, days=None,
                   start_date=None, end_date=None, fileprefix=''):
        """
        Initialize the fp_ohm class:
            1) initialize geopandas dataframe of concatenated hru_shapefiles
            2) initialize climate data using xarray
        :param ivars: list of vars to extract.  limited to ['tmax', 'tmin', 'ppt', 'rhmax', 'rhmin', 'ws', 'srad']
        :param iptpath: Input path, downloaded gridmet files will be saved here.  Must also contain shapefile used to
                        generate the weights file.
        :param optpath: Path to write newly generated netcdf file containing values for vars for each HRU
        :param weights_file: Weights file, based on shapefile in iptpath that was used to generate weights file
        :param etype: extraction time, can be days or date
        :param days: if extraction type == days, then # days to extract from most recently available date
        :param start_date: if extraction type date then start date in 'YYYY-MM-DD"
        :param end_date: if extraction type date then end date in 'YYYY-MM-DD"
        :param fileprefix: String to add to both downloaded gridment data and mapped hru file

        :return: success or failure
        """
        self.partial = partial
        for var in ivars:
            self.vars.append(var)
        self.iptpath = Path(iptpath)
        if self.iptpath.exists():
            print(f'input path exists {self.iptpath}', flush=True)
        else:
            sys.exit(f'Input Path does not exist: {self.iptpath} - EXITING')

        self.optpath = Path(optpath)
        if self.optpath.exists():
            print('output path exists', flush=True)
        else:
            sys.exit(f'Output Path does not exist: {self.optpath} - EXITING')

        self.wghts_file = Path(weights_file)
        if self.wghts_file.exists():
            print('weights file exists', self.wghts_file, flush=True)
        else:
            sys.exit(f'Weights file does not exist: {self.wghts_file} - EXITING')
        self.type = etype
        self.numdays = days
        self.start_date = start_date
        self.end_date = end_date
        self.fileprefix = fileprefix

        print(Path.cwd())
        if self.type == 'date':
            print(f'start_date: {self.start_date} and end_date: {self.end_date}', flush=True)
        else:
            print(f'number of days: {self.numdays}', flush=True)
        # glob.glob produces different results on Win and Linux. Adding sorted makes result consistent
        # filenames = sorted(glob.glob('*.shp'))
        # use pathlib glob
        # glob is here because original nhm had multiple shapefiles
        filenames = sorted(self.iptpath.glob('*.shp'))
        self.gdf = pd.concat([geopandas.read_file(f) for f in filenames], sort=True).pipe(geopandas.GeoDataFrame)
        self.gdf.reset_index(drop=True, inplace=True)

        print(f'the shapefile filenames read: {filenames}', flush=True)
        print(f'the shapefile header is: {self.gdf.head()}', flush=True)

        # self.num_hru = len(self.gdf.index)

        if self.type == 'date':
            self.numdays = ((self.end_date - self.start_date).days + 1)

        # Download netcdf subsetted data
        ncfile = []
        for var in self.vars:
            try:
                if var == 'tmax':
                    # Maximum Temperature
                    self.str_start, url, params = get_gm_url(self.type, 'tmax', self.numdays,
                                                             self.start_date, self.end_date)
                    self.write_extract_file(var, ncfile, url, params)
                    self.dstmax = xr.open_dataset(ncfile[-1])
                elif var == 'tmin':
                    # Minimum Temperature
                    self.str_start, url, params = get_gm_url(self.type, 'tmin', self.numdays,
                                                             self.start_date, self.end_date)
                    self.write_extract_file(var, ncfile, url, params)
                    self.dstmin = xr.open_dataset(ncfile[-1])
                elif var == 'ppt':
                    # Precipitation
                    self.str_start, url, params = get_gm_url(self.type, 'ppt', self.numdays,
                                                             self.start_date, self.end_date)
                    self.write_extract_file(var, ncfile, url, params)
                    self.dsppt = xr.open_dataset(ncfile[-1])
                elif var == 'rhmax':
                    # Maximum Relative Humidity
                    self.str_start, url, params = get_gm_url(self.type, 'rhmax', self.numdays,
                                                             self.start_date, self.end_date)
                    self.write_extract_file(var, ncfile, url, params)
                    self.dsrhmax = xr.open_dataset(ncfile[-1])
                elif var == 'rhmin':
                    # Minimum Relative Humidity
                    self.str_start, url, params = get_gm_url(self.type, 'rhmin', self.numdays,
                                                             self.start_date, self.end_date)
                    self.write_extract_file(var, ncfile, url, params)
                    self.dsrhmin = xr.open_dataset(ncfile[-1])
                elif var == 'ws':
                    # Mean daily Wind Speed
                    self.str_start, url, params = get_gm_url(self.type, 'ws', self.numdays,
                                                             self.start_date, self.end_date)
                    self.write_extract_file(var, ncfile, url, params)
                    self.dsws = xr.open_dataset(ncfile[-1])
                elif var == 'srad':
                    # Surface downwelling shortwave flux in air
                    self.str_start, url, params = get_gm_url(self.type, 'srad', self.numdays,
                                                             self.start_date, self.end_date)
                    self.write_extract_file(var, ncfile, url, params)
                    self.dssrad = xr.open_dataset(ncfile[-1])
            except HTTPError as http_err:
                print(f'HTTP error occured: {http_err}', flush=True)
                if self.numdays == 1:
                    sys.exit("numdays == 1: Gridmet not updated - EXITING")
                else:
                    sys.exit("GridMet not available or a bad request - EXITING")
            except Exception as err:
                sys.exit(f'Other error occured: {err}')
            else:
                print(f'Gridmet variable {var} retrieved: {ncfile[-1]}', flush=True)
        # self.ds = xr.open_mfdataset(ncfile, combine='by_coords')

        self.lat_h = self.dstmax['lat']
        self.lon_h = self.dstmax['lon']
        self.time_h = self.dstmax['day']

        ts = self.dstmax.sizes
        self.dayshape = ts['day']
        self.lonshape = ts['lon']
        self.latshape = ts['lat']

        # if self.type == 'days':
        print(f'Gridmet returned days = {self.dayshape} and expected number of days {self.numdays}', flush=True)
        if self.dayshape == self.numdays:
            return True
        else:
            print('returned and expected days not equal', flush=True)
            return False

    def run_weights(self):

        # read the weights file
        d_tmax = d_tmin = d_ppt = d_rhmax = d_rhmin = d_ws = d_srad = None
        d_flt_tmax = d_flt_tmin = d_flt_ppt = None
        d_flt_rhmin = d_flt_rhmax = d_flt_ws = d_flt_srad = None

        wght_uofi = pd.read_csv(self.wghts_file)
        # grab the hru_id from the weights file and use as identifier below
        self.wghts_id = wght_uofi.columns[1]

        # this geodataframe merges all hru-ids dissolves so the length of the index
        # equals the number of hrus
        self.gdf1 = self.gdf.sort_values(self.wghts_id).dissolve(by=self.wghts_id)
        self.num_hru = len(self.gdf1.index)

        # group by the weights_id for processing
        self.unique_hru_ids = wght_uofi.groupby(self.wghts_id)

        print('finished reading weight file', flush=True)

        # intialize numpy arrays to store climate vars
        self.np_tmax = zeros((self.numdays, self.num_hru))
        self.np_tmin = zeros((self.numdays, self.num_hru))
        self.np_ppt = zeros((self.numdays, self.num_hru))
        self.np_rhmax = zeros((self.numdays, self.num_hru))
        self.np_rhmin = zeros((self.numdays, self.num_hru))
        self.np_ws = zeros((self.numdays, self.num_hru))
        self.np_srad = zeros((self.numdays, self.num_hru))

        tindex = asarray(self.gdf1.index, dtype=np.int)

        for day in arange(self.numdays):
            print(f'Processing day: {day}', flush=True)

            d_tmax = zeros(self.num_hru)
            d_flt_tmax = self.dstmax[self.gmss_vars['tmax']].values[day, :, :].flatten(order='K')

            d_tmin = zeros(self.num_hru)
            d_flt_tmin = self.dstmin[self.gmss_vars['tmin']].values[day, :, :].flatten(order='K')

            d_ppt = zeros(self.num_hru)
            d_flt_ppt = self.dsppt[self.gmss_vars['ppt']].values[day, :, :].flatten(order='K')

            d_rhmax = zeros(self.num_hru)
            d_flt_rhmax = self.dsrhmax[self.gmss_vars['rhmax']].values[day, :, :].flatten(order='K')

            d_rhmin = zeros(self.num_hru)
            d_flt_rhmin = self.dsrhmin[self.gmss_vars['rhmin']].values[day, :, :].flatten(order='K')

            d_ws = zeros(self.num_hru)
            d_flt_ws = self.dsws[self.gmss_vars['ws']].values[day, :, :].flatten(order='K')

            if 'srad' in self.vars:
                d_srad = zeros(self.num_hru)
                d_flt_srad = self.dssrad[self.gmss_vars['srad']].values[day, :, :].flatten(order='K')

            for i in arange(len(tindex)):
                if not self.partial:
                    weight_id_rows = self.unique_hru_ids.get_group(tindex[i])
                    tw = weight_id_rows.w.values
                    tgid = weight_id_rows.grid_ids.values
                    # tmask, tgid, tw = getweights(index[i], gid, hid, w)

                    d_tmax[i] = getaverage(d_flt_tmax[tgid] - 273.15, tw)
                    d_tmin[i] = getaverage(d_flt_tmin[tgid] - 273.15, tw)
                    d_ppt[i] = getaverage(d_flt_ppt[tgid], tw)
                    d_rhmax[i] = getaverage(d_flt_rhmax[tgid], tw)
                    d_rhmin[i] = getaverage(d_flt_rhmin[tgid], tw)
                    d_ws[i] = getaverage(d_flt_ws[tgid], tw)

                    if 'srad' in self.vars:
                        d_srad[i] = getaverage(d_flt_srad[tgid], tw)
                else:
                    try:
                        weight_id_rows = self.unique_hru_ids.get_group(tindex[i])
                        tw = weight_id_rows.w.values
                        tgid = weight_id_rows.grid_ids.values
                        # if one var is nan all are nan. getaverage returns nan if 1 value is nan
                        # np_get_wval return nan if all values are nan otherwise returns the
                        # weighted masked val.  So assumption here is return a value for partially
                        # weighted HRUs
                        nanvar = False
                        if np.isnan(getaverage(d_flt_tmax[tgid], tw)):
                            nanvar = True

                        if nanvar:
                            d_tmax[i] = np_get_wval(d_flt_tmax[tgid] - 273.15, tw, tindex[i])
                            d_tmin[i] = np_get_wval(d_flt_tmin[tgid] - 273.15, tw, tindex[i])
                            d_ppt[i] = np_get_wval(d_flt_ppt[tgid], tw, tindex[i])
                            d_rhmax[i] = np_get_wval(d_flt_rhmax[tgid], tw, tindex[i])
                            d_rhmin[i] = np_get_wval(d_flt_rhmin[tgid], tw, tindex[i])
                            d_ws[i] = np_get_wval(d_flt_ws[tgid], tw, tindex[i])

                        else:
                            d_tmax[i] = getaverage(d_flt_tmax[tgid] - 273.15, tw)
                            d_tmin[i] = getaverage(d_flt_tmin[tgid] - 273.15, tw)
                            d_ppt[i] = getaverage(d_flt_ppt[tgid], tw)
                            d_rhmax[i] = getaverage(d_flt_rhmax[tgid], tw)
                            d_rhmin[i] = getaverage(d_flt_rhmin[tgid], tw)
                            d_ws[i] = getaverage(d_flt_ws[tgid], tw)

                        if 'srad' in self.vars:
                            if nanvar:
                                d_srad[i] = np_get_wval(d_flt_srad[tgid], tgid, tw, tindex[i])
                            else:
                                d_srad[i] = getaverage(d_flt_srad[tgid], tw)
                    except KeyError:
                        # This except block protects against HRUs that are completely
                        # outside the footprint of Gridmet.  If so, they will have no value
                        # in the weights file and so return default value.
                        d_tmax[i] = netCDF4.default_fillvals['f8']
                        d_tmin[i] = netCDF4.default_fillvals['f8']
                        d_ppt[i] = netCDF4.default_fillvals['f8']
                        d_rhmax[i] = netCDF4.default_fillvals['f8']
                        d_rhmin[i] = netCDF4.default_fillvals['f8']
                        d_ws[i] = netCDF4.default_fillvals['f8']
                        if 'srad' in self.vars:
                            d_srad[i] = netCDF4.default_fillvals['f8']

                    if i % 10000 == 0:
                        print(f'    Processing hru {i}', flush=True)


            self.np_tmax[day, :] = d_tmax[:]
            self.np_tmin[day, :] = d_tmin[:]
            self.np_ppt[day, :] = d_ppt[:]
            self.np_rhmax[day, :] = d_rhmax[:]
            self.np_rhmin[day, :] = d_rhmin[:]
            self.np_ws[day, :] = d_ws[:]
            if 'srad' in self.vars:
                self.np_srad[day, :] = d_srad[:]

        # self.dstmax.close()

    def finalize(self):
        print(Path.cwd(), flush=True)

        ncfile = Dataset(
            self.optpath / (self.fileprefix + 'climate_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc'),
            mode='w', format='NETCDF4_CLASSIC')

        def getxy(pt):
            return pt.x, pt.y

        centroidseries = self.gdf1.geometry.centroid
        tlon, tlat = [list(t) for t in zip(*map(getxy, centroidseries))]

        # Global Attributes
        ncfile.Conventions = 'CF-1.8'
        ncfile.featureType = 'timeSeries'
        ncfile.history = ''

        # Create dimensions
        for dim in ncfile.dimensions.items():
            print(dim, flush=True)

        # Create Variables
        time = ncfile.createVariable('time', 'f4', ('time',))
        time.long_name = 'time'
        time.standard_name = 'time'
        time.units = 'days since ' + self.str_start
        time.calendar = 'standard'
        time[:] = arange(0, self.numdays, dtype=np.float)

        hru = ncfile.createVariable('hruid', 'i', ('hruid',))
        hru.cf_role = 'timeseries_id'
        hru.long_name = 'local model hru id'
        hru[:] = np.asarray(self.gdf1.index)

        lat = ncfile.createVariable('hru_lat', dtype(float32).char, ('hruid',))
        lat.long_name = 'Latitude of HRU centroid'
        lat.units = 'degrees_north'
        lat.standard_name = 'hru_latitude'
        lat[:] = tlat

        lon = ncfile.createVariable('hru_lon', dtype(float32).char, ('hruid',))
        lon.long_name = 'Longitude of HRU centroid'
        lon.units = 'degrees_east'
        lon.standard_name = 'hru_longitude'
        lon[:] = tlon

        for var in self.vars:
            if var == 'tmax':
                tmax = ncfile.createVariable('tmax', dtype(float32).char, ('time', 'hruid'))
                tmax.long_name = 'Maximum daily air temperature'
                tmax.units = 'degree_Celsius'
                tmax.standard_name = 'maximum_daily_air_temperature'
                tmax.fill_value = default_fillvals['f8']
                tmax[:, :] = self.np_tmax[:, :]

            elif var == 'tmin':
                tmin = ncfile.createVariable('tmin', dtype(float32).char, ('time', 'hruid'))
                tmin.long_name = 'Minimum daily air temperature'
                tmin.units = 'degree_Celsius'
                tmin.standard_name = 'minimum_daily_air_temperature'
                tmin.fill_value = default_fillvals['f8']
                tmin[:, :] = self.np_tmin[:, :]

            elif var == 'ppt':
                prcp = ncfile.createVariable('prcp', dtype(float32).char, ('time', 'hruid'))
                prcp.long_name = 'Daily Accumulated Precipitation'
                prcp.units = 'mm'
                prcp.standard_name = 'prcp'
                prcp.fill_value = default_fillvals['f8']
                prcp[:, :] = self.np_ppt[:, :]

            elif var == 'rhmax':
                rhmax = ncfile.createVariable('rhmax', dtype(float32).char, ('time', 'hruid'))
                rhmax.long_name = 'Daily Maximum Relative Humidity'
                rhmax.units = 'percent'
                rhmax.standard_name = 'rhmax'
                rhmax.fill_value = default_fillvals['f8']
                rhmax[:, :] = self.np_rhmax[:, :]

            elif var == 'rhmin':
                rhmin = ncfile.createVariable('rhmin', dtype(float32).char, ('time', 'hruid'))
                rhmin.long_name = 'Daily Maximum Relative Humidity'
                rhmin.units = 'percent'
                rhmin.standard_name = 'rhmin'
                rhmin.fill_value = default_fillvals['f8']
                rhmin[:, :] = self.np_rhmin[:, :]

            elif var == 'ws':
                ws = ncfile.createVariable('ws', dtype(float32).char, ('time', 'hruid'))
                ws.long_name = 'Daily Mean Wind Speed'
                ws.units = 'm/s'
                ws.standard_name = 'ws'
                ws.fill_value = default_fillvals['f8']
                ws[:, :] = self.np_ws[:, :]

            elif var == 'srad':
                srad = ncfile.createVariable('srad', dtype(float32).char, ('time', 'hruid'))
                srad.long_name = 'surface_downwelling_shortwave_flux_in_air'
                srad.units = 'W m-2'
                srad.standard_name = 'srad'
                srad.fill_value = default_fillvals['f8']
                srad[:, :] = self.np_srad[:, :]

        hum = ncfile.createVariable('humidity', dtype(float32).char, ('time', 'hruid'))
        hum.long_name = 'Daily mean relative humidity'
        hum.units = 'percent'
        hum.standard_name = 'rhavg'
        hum.fill_value = default_fillvals['f8']
        hum[:, :] = (self.np_rhmax[:, :] + self.np_rhmin[:, :])/2.0

        ncfile.close()
        print("dataset is closed", flush=True)

    def setnumdays(self, num_d):
        self.numdays = num_d
