from argparse import ArgumentParser, Namespace
from cartopy import geodesic
import cdsapi
from Config import Config
from datetime import datetime, timedelta, timezone
import logging
from os import path
import netCDF4 as nc
import numpy
from shapely import Point

logger = logging.getLogger(__name__)

class ERA5:
    '''
    Common properties of ERA5 grids
    '''
    def __init__(self):

        self.lat_index = -1
        self.lat_len = 0
        self.lat_min = 90.
        self.lat_max  = -90.

        self.lon_index = -1
        self.lon_len = 0
        self.lon_min = 180.
        self.lon_max = -180.

        self.time_index = -1
        self.time_min = 0
        self.time_len = 0
        self.time_max = 0

        self.time: numpy.array = None
        self.lats: numpy.array = None
        self.lons: numpy.array = None
        self.variables: dict[str,numpy.array] = {}

    def calculate_ERA_GridTime(observation_time: datetime) -> float:
        '''
        Find ERA 5 grid time nearest to time of observation

        Parameters
        ----------
        observation_time: datetime
            time of observation

        Returns
        -------
        era_grid_timestamp: float
            ERA5 grid timestamp closest to time
        '''
        # level = len(inspect.stack())
        # print("{}ERA5: calculate_ERA_GridTime".format("\t" * level))

        dtout = datetime(observation_time.year, observation_time.month, observation_time.day, observation_time.hour , tzinfo=timezone.utc)

        if (30 <= observation_time.minute):
            dtout += timedelta(hours=1)

        # print("{}Pixel time {}".format("\t" * (level + 1), pixeltime.strftime("%Y%m%d %H:%M:%S")))
        # print("{}Nearest ERA5 time {}".format("\t" * (level + 1), dtout.strftime("%Y%m%d %H:%M:%S")))

        return dtout.timestamp()

    def calculate_ERA5_NearestGridPoint(observation_point: Point) -> Point:
        '''
        Find ERA 5 grid point nearest to observation point

        Parameters
        ----------
        observation_point: Point
            Location of observation 

        Returns
        -------
        era_grid_point: Point
            Nearest ERA5 grid point
        '''
        # level = len(inspect.stack())
        # print("{}ERA5: calculate_ERA5_NearestGridPoint".format("\t" * level))

        g = geodesic.Geodesic()
        lon = int(observation_point.x)
        lat = int(observation_point.y)
        if lon < 0: lon -= 1
        if lat < 0: lat -= 1
        nlat = None
        nlon = None
        mindistance = 1000000
        for n in range(4):
            for t in range(4):
                dx = lon + n * 0.25
                dy = lat + t * 0.25
                [[distance, _, _]] = g.inverse([observation_point.x, observation_point.y], [dx,  dy])
                if (distance < mindistance):
                    nlat = dy
                    nlon = dx
                    mindistance = distance

        # print("{}Mine: Lon: {} Lat: {}".format("\t" * (level + 1), observation_point.x, observation_point.y))
        # print("{}Nearast ERA5 grid point: Lon: {} Lat: {}".format("\t" * (level + 1), nlon, nlat))
        return Point(nlon,nlat)

    def convert_lat_to_index(self, y: float) -> None:
        '''
        Convert grid latitude to grid index. On exit value self.lat_index is set to index or -1 if not found

        Parameters
        ----------
        y: float
            Latitude
        '''
        self.lat_index = -1
        if (y < self.lat_min): return
        if (self.lat_max < y): return
        for i in range(self.lat_len):
            if self.lats[i] == y:
                self.lat_index = i
                return

    def convert_lon_to_index(self, x: float) -> None:
        '''
        Convert grid longitude to grid index. On exit value self.lon_index is set to index or -1 if not found

        Parameters
        ----------
        x: float
            Longitude
        '''
        self.lon_index = -1
        if (x < self.lon_min): return
        if (self.lon_max < x): return
        for i in range(self.lon_len):
            if self.lons[i] == x:
                self.lon_index = i
                return

    def convert_timestamp_to_index(self, t: float) -> None:
        '''
        Convert grid time to grid index. On exit value self.time_index is set to index or -1 if not found

        Parameters
        ----------
        t: float
            Timestamp
        '''
        self.time_index = -1
        if (t < self.time_min): return
        if (self.time_max < t): return
        for i in range(self.time_len):
            if self.time[i] == t:
                self.time_index = i
                return

class ERA5_SingleLevel(ERA5):
    '''
    Class used to manipulate ERA5 Reanalysis Single Level Grid.
    '''
    def __init__(self):
        super().__init__()

    @classmethod
    def download(cls, request: dict, target: str):
        '''
        Factory method downloading data into a file and returning populated object

        Parameters
        ----------
        request: dict
            Content of the ERA5 dataset in the format of cdsapi
        target: str
            Destination of the downloaded data set

        Returns
        -------
        ERA5_SingleLevel object populated with data
        '''
        dataset = "reanalysis-era5-single-levels"
        client = cdsapi.Client()
        client.retrieve(dataset, request, target)
        cls = ERA5_SingleLevel.from_netCDF(target)
        return cls

    @classmethod
    def download_ERA5_hourly_BLH(cls, source: str, x: float, y:float):
        '''
        Download ERA5 Boundary Layer Height from 2018-01-01 00Z to 2019:12:31 23Z
        '''
        path_ERA5_SL = path.join(Config.ERA5_folder, source, "ERA5_BLH_2018_2019_{}.nc".format(source))
        if not path.exists(path_ERA5_SL):
            request = {
                    "product_type": ["reanalysis"],
                    "variable": ["boundary_layer_height"],
                    "year": ["2018","2019"],
                    "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12" ],
                    "day": ["01", "02", "03","04", "05", "06","07", "08", "09", "10", 
                            "11", "12","13", "14", "15","16", "17", "18", "19", "20",
                            "21", "22", "23", "24","25", "26", "27","28", "29", "30",
                            "31"],
                    "time": [
                        "00:00", "01:00", "02:00","03:00", "04:00", "05:00","06:00", "07:00", "08:00","09:00", "10:00", "11:00",
                        "12:00", "13:00", "14:00","15:00", "16:00", "17:00","18:00", "19:00", "20:00","21:00", "22:00", "23:00"
                    ],
                    "data_format": "netcdf",
                    "download_format": "unarchived",
                    "area": [y, x, y, x]
                }
            ERA5_SingleLevel.download(request, path_ERA5_SL)

    @classmethod
    def from_netCDF(cls, filename: str):
        '''
        Factory method loading dataset from ERA5 reanalysis netcdf file and update state of the class

        Parameters
        ----------
        filename: str
            Path to the netcdf file

        Returns
        -------
        obj: ERA5_SingleLevel
            Object populated with data from netCDF file
        '''
        if (not(logger is None)): logger.info("ERA5_SingleLevel: from_netCDF File: {}".format(filename))
        ds  = nc.Dataset(filename)
        cls = ERA5_SingleLevel()
        for k,v in ds.variables.items():
            if (not(logger is None)): logger.info("Available variable: {}".format(k))
            if k == "number" : continue
            elif k == "expver" : continue
            elif k == "valid_time" : cls.time = v[:]
            elif k == "latitude" : cls.lats = v[:]
            elif k == "longitude": cls.lons = v[:]
            else: cls.variables[k] = v[:, :, :]
        ds.close()

        cls.lat_min = numpy.min(cls.lats)
        cls.lat_max = numpy.max(cls.lats)
        cls.lat_len = len(cls.lats)
        cls.lon_min = numpy.min(cls.lons)
        cls.lon_max = numpy.max(cls.lons)
        cls.lon_len = len(cls.lons)
        cls.time_min = numpy.min(cls.time)
        cls.time_max = numpy.max(cls.time)
        cls.time_len = len(cls.time)
        if (not(logger is None)): logger.info("t:[{} - {}] lon:[{} - {}] lat:[{} - {}]".format(cls.time_min, cls.time_max, cls.lon_min, cls.lon_max, cls.lat_min, cls.lat_max))

        return cls

    def get(self, variablename: str, timeindex: int, latindex: int, lonindex: int) -> float:
        '''
        Get a variable at grid point.
        Known variables are: blh, mslp, sp, z

        Parameters
        ----------
        variablename: str
            Known variables are: blh, mslp, sp, z, u10, v10
        timeindex: int
            Grid index of time
        latindex: int
            Grid index of latitude
        lonindex: int
            Grid index of longitude

        Returns
        -------
        Value of the variable or exception if parameters outside of this dataset
        '''
        if (not(variablename in self.variables)): raise Exception("Variable {} does not exist".format(variablename))

        if (timeindex < 0): raise Exception("Time index {} cannot be negative".format(timeindex))
        if (len(self.time) < timeindex): raise Exception("Time index {} outside of grid {}".format(timeindex, len(self.time)))

        if (latindex < 0): raise Exception("Latitude index {} cannot be negative".format(latindex))
        if (len(self.lats) < latindex): raise Exception("Latitude index {} outside of grid {}".format(latindex, len(self.lats)))

        if (lonindex < 0): raise Exception("Longitude index {} cannot be negative".format(lonindex))
        if (len(self.lons) < lonindex): raise Exception("Longitude index {} outside of grid {}".format(lonindex, len(self.lons)))

        return self.variables[variablename][timeindex, latindex, lonindex]

    @staticmethod
    def get_filepath(algorithm: str, source: str, orbit:str) -> str:
        return path.join(Config.ERA5_folder, source, "ERA5_{}_{}_SingleLevel.nc".format(algorithm, orbit))


class ERA5_PressureLevels(ERA5):
    '''
    Class used to manipulate ERA5 Reanalysis Pressure Level Grids.    
    https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=overview
    DOI: 10.24381/cds.bd0915c6
    '''

    def __init__(self):
        super().__init__()
        self.p_index = -1
        self.p_len = 0
        self.p_min = None
        self.p_max = None
        self.pressure: numpy.array = None

    def calculate_pressure_index(self, p: float) -> None:
        self.p_index = -1
        for i in range(self.p_len):
            if self.pressure[i] == p:
                self.p_index = i
                break

    @classmethod
    def download(cls, request: dict, target: str):
        '''
        Download ERA5 reanalysis on pressure levels as specified in the dictionary and write to the target file

        Parameters
        ----------
        request: dict
            Content of the ERA5 dataset in the format of cdsapi
        target: str
            Destination of the downloaded data set

        Returns
        -------
        None
        '''
        dataset = "reanalysis-era5-pressure-levels"
        client = cdsapi.Client()
        client.retrieve(dataset, request, target)
        cls = ERA5_PressureLevels.from_netCDF(target)
        return cls

    @classmethod
    def from_netCDF(cls, filename: str):
        '''
        Reconstitute this data set from netCDF4 file

        Parameters
        ----------
        filename: str
            Path to the netcdf file

        '''
        if (not(logger is None)): logger.info("ERA5_PressureLevels: from_netCDF: File: {}".format(filename))
        ds  = nc.Dataset(filename)
        cls = ERA5_PressureLevels()
        cls.variables = {}
        for k,v in ds.variables.items():
            if (not(logger is None)): logger.info("Available variable: {}".format(k))
            if k == "number" : continue
            elif k == "expver" : continue
            elif k == "valid_time" : cls.time = v[:]
            elif k == "latitude" : cls.lats = v[:]
            elif k == "longitude": cls.lons = v[:]
            elif k == "pressure_level": 
                cls.pressure = v[:]
                cls.p_len = len(v)
            else: cls.variables[k] = v[:, :, :, :]
        ds.close()

        cls.lat_index = -1
        cls.lat_len = len(cls.lats)
        cls.lat_min = numpy.min(cls.lats)
        cls.lat_max = numpy.max(cls.lats)

        cls.lon_index = -1
        cls.lon_len = len(cls.lons)
        cls.lon_min = numpy.min(cls.lons)
        cls.lon_max = numpy.max(cls.lons)

        cls.p_index = -1
        cls.p_len = len(cls.pressure)
        cls.p_min = numpy.min(cls.pressure)
        cls.p_max = numpy.max(cls.pressure)
        
        cls.time_index = -1
        cls.time_len = len(cls.time)
        cls.time_min = numpy.min(cls.time)
        cls.time_max = numpy.max(cls.time)

        if (not(logger is None)): logger.info("t:[{} - {}] p:[{} - {}]  lon:[{} - {}] lat:[{} - {}]".format(cls.time_min, cls.time_max, cls.p_max, cls.p_min, cls.lon_min, cls.lon_max, cls.lat_min, cls.lat_max))

        return cls

    def get(self, variablename: str, timeindex: int, pressureindex: int, latindex: int, lonindex: int) -> float:
        '''
        Get a variable at grid point.
        Known variables are: p, u, v, z

        Parameters
        ----------
        variablename: str
            Known variables are: blh, mslp, sp, z
        timeindex: int
            Grid index of time
        pressureindex: int
            Grid index of pressure
        latindex: int
            Grid index of latitude
        lonindex: int
            Grid index of longitude

        Returns
        -------
        Value of the variable or exception if parameters outside of this dataset
        '''
        if (not(variablename in self.variables)): raise Exception("Variable {} does not exist".format(variablename))

        if (timeindex < 0): raise Exception("Time index {} cannot be negative".format(timeindex))
        if (len(self.time) < timeindex): raise Exception("Time index {} outside of grid {}".format(timeindex, len(self.time)))

        if (pressureindex < 0): raise Exception("Pressure index {} cannot be negative".format(pressureindex))
        if (len(self.pressure) < pressureindex): raise Exception("Pressure index {} outside of grid {}".format(pressureindex, len(self.pressure)))

        if (latindex < 0): raise Exception("Latitude index {} cannot be negative".format(latindex))
        if (len(self.lats) < latindex): raise Exception("Latitude index {} outside of grid {}".format(latindex, len(self.lats)))

        if (lonindex < 0): raise Exception("Longitude index {} cannot be negative".format(lonindex))
        if (len(self.lons) < lonindex): raise Exception("Longitude index {} outside of grid {}".format(lonindex, len(self.lons)))

        return self.variables[variablename][timeindex, pressureindex, latindex, lonindex]

    @staticmethod
    def get_filepath(algorithm: str, source: str, orbit:str) -> str:
        return path.join(Config.ERA5_folder, source, "ERA5_{}_{}_PressureLevels.nc".format(algorithm, orbit))

def _handler_download(args: Namespace):
    '''
    Handler for download of ERA5 data on pressure levels
    '''
    from Config import Config
    from Source import create_source
    from TROPOMI import TROPOMI, TROPOMI_for_orbit

    config = Config()
    source = create_source(args.source)
    Config.create_directory_structure_ERA5(source.case_name)

    if args.algorithm == "CSF" :
        Config.create_log("ERA5_CSF_{}.log".format(args.source))

        path_ERA5_PL = ERA5_PressureLevels.get_filepath("CSF",source.case_name, args.orbit)
        path_ERA5_SL = ERA5_SingleLevel.get_filepath("CSF", source.case_name, args.orbit)

        if path.exists(path_ERA5_PL) and path.exists(path_ERA5_SL): return

        e = False
        for p in Config.TROPOMI_Processors:
            if TROPOMI.exists_locally(config, args.orbit, p):
                e = True
                break

        if not e:
            message = "Unable to find required TROPOMI file for orbit {} . Please download this file first.".format(args.orbit)
            if not logger is None: logger.error(message)
            print(message)
            return

        tropomi = TROPOMI_for_orbit(args.orbit, args.processor)
        [_, _, dt, _] = tropomi.get_pixel_for_source(source.xy)
        grid_time = ERA5.calculate_ERA_GridTime(dt)
        
        dt1 = datetime.fromtimestamp(grid_time, tz=timezone.utc)
        str_year = "{:04d}".format(dt1.year)
        str_month = "{:02d}".format(dt1.month)
        str_day = "{:02d}".format(dt1.day)
        str_hour = "{:02d}:00".format(dt1.hour)

        if not path.exists(path_ERA5_PL):
            request = {
                "product_type": ["reanalysis"],
                "variable": [
                    "geopotential",
                    "u_component_of_wind",
                    "v_component_of_wind"
                ],
                "year": [str_year],
                "month": [str_month],
                "day": [str_day],
                "time": [str_hour],
                "pressure_level": [
                    "1", "2", "3",
                    "5", "7", "10",
                    "20", "30", "50",
                    "70", "100", "125",
                    "150", "175", "200",
                    "225", "250", "300",
                    "350", "400", "450",
                    "500", "550", "600",
                    "650", "700", "750",
                    "775", "800", "825",
                    "850", "875", "900",
                    "925", "950", "975",
                    "1000"
                ],
                "data_format": "netcdf",
                "download_format": "unarchived",
                "area": [source.gridpoint.y, source.gridpoint.x, source.gridpoint.y, source.gridpoint.x]        
            }
            _ = ERA5_PressureLevels.download(request,  path_ERA5_PL)

        if not path.exists(path_ERA5_SL):
            request = {
                "product_type": ["reanalysis"],
                "variable": ["boundary_layer_height", "surface_pressure"],
                "year": [str_year],
                "month": [str_month],
                "day": [str_day],
                "time": [str_hour],
                "data_format": "netcdf",
                "download_format": "unarchived",
                "area": [source.gridpoint.y, source.gridpoint.x, source.gridpoint.y, source.gridpoint.x]
            }
            _ = ERA5_SingleLevel.download(request, path_ERA5_SL)

def _handler_list(args: Namespace):
    '''
    Handle for listing content of files
    '''
    from Config import Config
    from Constants import GRAVITY
    from Source import create_source
    from TROPOMI import TROPOMI_for_orbit
    Config.create_log("ERA5_CSF_{}_{}.log".format(args.source, args.orbit))
    source = create_source(args.source)
    orbit = args.orbit
    tropomi = TROPOMI_for_orbit(args.orbit, args.processor)

    if tropomi is None:
        m = "Error: Missing TROPOMI data for orbit: {}. Please download first.".format(args.orbit)
        print(m)
        if not logger is None: logger.error(m)
        return

    [scan, pixel, _, dt] = tropomi.get_pixel_for_source(source.xy)
    if dt is None:
        m = "Error: TROPOMI orbit: {} does not contain source {}.".format(args.orbit, args.source)
        print(m)
        if not logger is None: logger.error(m)
        return

    path_ERA5_SL = ERA5_SingleLevel.get_filepath("CSF", source.case_name, orbit)
    
    if (not(path.exists(path_ERA5_SL))):
        m = "Error: This function expects to use file {}. Pleas make sure that ERA5 file is downloaded first.".format(path_ERA5_SL)
        print(m)
        if not logger is None: logger.error(m)
        return
    
    era5_sl = ERA5_SingleLevel.from_netCDF(path_ERA5_SL)
    
    path_ERA5_PL = ERA5_PressureLevels.get_filepath("CSF", source.case_name, orbit)
    if (not(path.exists(path_ERA5_PL))):
        m = "Error: This function expects to use file {}. Pleas make sure that ERA5 file is downloaded first.".format(path_ERA5_PL)
        print(m)
        if not logger is None: logger.error(m)
        return
    
    era5_p = ERA5_PressureLevels.from_netCDF(path_ERA5_PL)

    grid_time = ERA5.calculate_ERA_GridTime(dt)
    era5_p.convert_timestamp_to_index(grid_time)
    era5_p.convert_lat_to_index(source.gridpoint.y)
    era5_p.convert_lon_to_index(source.gridpoint.x)

    era5_sl.convert_timestamp_to_index(grid_time)
    era5_sl.convert_lat_to_index(source.gridpoint.y)
    era5_sl.convert_lon_to_index(source.gridpoint.x)

    blh = era5_sl.get("blh", era5_sl.time_index, era5_sl.lat_index, era5_sl.lon_index) 
    sp = era5_sl.get("sp", era5_sl.time_index, era5_sl.lat_index, era5_sl.lon_index) / 100.0

    p = era5_p.pressure
    z = era5_p.variables["z"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]
    u = era5_p.variables["u"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]
    v = era5_p.variables["v"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]

    messages = []
    messages.append("TROPOMI orbit: {}".format(args.orbit))
    messages.append("\tscan: {} pixel: {} nearest hour: {}Z".format(scan, pixel, dt.strftime("%Y-%m-%d:%H")))
    messages.append("Single level variables")
    messages.append("\tBoundary Layer Height: {} (m), Surface Pressure: {} (hPa)".format(blh, sp))
    messages.append("Pressure level variables")
    for i in range(len(p)):
        messages.append("\tPressure: {} (hPa), Height: {} (m), u: {} (m/s), v: {} (m/s)".format(p[i], z[i] / GRAVITY, u[i], v[i]))

    print("\n".join(messages))
    if not logger is None:
        for m in messages:
            logger.info(m)    
    pass

def _handler_wind(args: Namespace) -> None:
    '''
    Test for orbit Hail Creek 09956
    '''
    from Config import Config
    from math import sqrt
    from Meteorology import Meteorology
    from Source import create_source
    from TROPOMI import TROPOMI_for_orbit
    Config.create_log("ERA5_CSF_{}_{}.log".format(args.source, args.orbit))
    source = create_source(args.source)
    orbit = args.orbit
    tropomi = TROPOMI_for_orbit(args.orbit, args.processor)

    if tropomi is None:
        m = "Error: Missing TROPOMI data for orbit: {}. Please download first.".format(args.orbit)
        print(m)
        if not logger is None: logger.error(m)
        return

    [_, _, _, dt] = tropomi.get_pixel_for_source(source.xy)
    if dt is None:
        m = "Error: TROPOMI orbit: {} does not contain source {}.".format(args.orbit, args.source)
        print(m)
        if not logger is None: logger.error(m)
        return

    path_ERA5_SL = ERA5_SingleLevel.get_filepath("CSF", source.case_name, orbit)
    
    if (not(path.exists(path_ERA5_SL))):
        m = "Error: This function expects to use file {}. Pleas make sure that ERA5 file is downloaded first.".format(path_ERA5_SL)
        print(m)
        if not logger is None: logger.error(m)
        return
    
    era5_sl = ERA5_SingleLevel.from_netCDF(path_ERA5_SL)
    
    path_ERA5_PL = ERA5_PressureLevels.get_filepath("CSF", source.case_name, orbit)
    if (not(path.exists(path_ERA5_PL))):
        m = "Error: This function expects to use file {}. Pleas make sure that ERA5 file is downloaded first.".format(path_ERA5_PL)
        print(m)
        if not logger is None: logger.error(m)
        return
    
    era5_p = ERA5_PressureLevels.from_netCDF(path_ERA5_PL)

    grid_time = ERA5.calculate_ERA_GridTime(dt)
    era5_p.convert_timestamp_to_index(grid_time)
    era5_p.convert_lat_to_index(source.gridpoint.y)
    era5_p.convert_lon_to_index(source.gridpoint.x)

    era5_sl.convert_timestamp_to_index(grid_time)
    era5_sl.convert_lat_to_index(source.gridpoint.y)
    era5_sl.convert_lon_to_index(source.gridpoint.x)

    blh = era5_sl.get("blh", era5_sl.time_index, era5_sl.lat_index, era5_sl.lon_index) 
    sp = era5_sl.get("sp", era5_sl.time_index, era5_sl.lat_index, era5_sl.lon_index) / 100.0

    p = era5_p.pressure
    z = era5_p.variables["z"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]
    u = era5_p.variables["u"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]
    v = era5_p.variables["v"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]

    wind = Meteorology.calculate_ERA5_PressureAveragedWind(p, z, u, v, blh, sp)
    messages = []
    messages.append("BLH pressure averaged wind. u: {} (m/s), v: {} (m/s), speed: {} (m/s)".format(wind[0], wind[1], sqrt(wind[0] * wind[0] + wind[1] * wind[1])))
    print("\n".join(messages))
    if not logger is None:
        for m in messages:
            logger.info(m)

if __name__ == '__main__':
    '''
    Entry point into this module
    '''
    from Source import Source
    parser = ArgumentParser(prog="ERA5", description="Utility to download ERA5 data for specific algorithm, source and orbit")
    subparsers = parser.add_subparsers()

    parser_download = subparsers.add_parser("download", help="Download data for algorithm, location and orbit")
    parser_download.add_argument("algorithm", type=str, choices=["CSF"], help="Selection of algorithm for which ERA files are downloaded")
    parser_download.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_download.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 04579")
    parser_download.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")    
    parser_download.set_defaults(func=_handler_download)

    parser_list = subparsers.add_parser("list", help="List data for algorithm, location and orbit")
    parser_list.add_argument("algorithm", type=str, choices=["CSF"], help="Selection of algorithm for which ERA files are downloaded")
    parser_list.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_list.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 04579")
    parser_list.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_list.set_defaults(func=_handler_list)

    parser_wind = subparsers.add_parser("wind", help="Calculate pressure averaged boundary layer wind for algorithm, location and orbit")
    parser_wind.add_argument("algorithm", type=str, choices=["CSF"], help="Selection of algorithm for which ERA files are downloaded")
    parser_wind.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_wind.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 04579")
    parser_wind.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_wind.set_defaults(func=_handler_wind)

    args = parser.parse_args()    

    if ("orbit" in args):
        if not(len(args.orbit) == 5) or not(args.orbit.isdigit()):
            print ("Orbit must be a five digit number: {}".format(args.orbit))
            exit()
    
    args.func(args)