from Config import Config
from datetime import datetime, timezone
from ERA5 import ERA5, ERA5_PressureLevels, ERA5_SingleLevel
import logging
from Meteorology import Meteorology
from os import path
from Source import create_source
import unittest

logger = logging.getLogger(__name__)

class Meteorology_test(unittest.TestCase):
    def __init__(self, methodName = "runTest"):
        super().__init__(methodName)

    @classmethod
    def setUpClass(cls):
        Config.create_log("Meteorology_test.log")

    @classmethod
    def tearDownClass(cls):
        pass

    def test_calculate_wind(self):
        '''
        Test for orbit Hail Creek 09956
        '''
        source = create_source("HailCreek")
        orbit = "09956"
        dt = datetime(2019, 9, 15, 4, tzinfo=timezone.utc)

        path_ERA5_SL = ERA5_SingleLevel.get_filepath("CSF", source.case_name, orbit)
        
        if (not(path.exists(path_ERA5_SL))):
            self.fail("Error: This test expects to use file {}. Pleas make sure that ERA5 test passes first.".format(path_ERA5_SL))
            return
        
        era5_sl = ERA5_SingleLevel.from_netCDF(path_ERA5_SL)
        
        path_ERA5_PL = ERA5_PressureLevels.get_filepath("CSF",source.case_name, orbit)
        if (not(path.exists(path_ERA5_PL))):
            self.fail("Error: This test expects to use file {}. Pleas make sure that ERA5 test passes first.".format(path_ERA5_PL))
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

        expected_u = -5.6952834868095294
        expected_v = -1.0697813437018595
        self.assertEqual(expected_u, wind[0])
        self.assertEqual(expected_v, wind[1])

    def test_calculate_azimuth_meteorology(self):
        self.assertEqual(0, Meteorology.calculate_azimuth_degree_meteorology(0, -1))
        self.assertEqual(90, Meteorology.calculate_azimuth_degree_meteorology(-1, 0))
        self.assertEqual(180, Meteorology.calculate_azimuth_degree_meteorology(0, 1))
        self.assertEqual(270, Meteorology.calculate_azimuth_degree_meteorology(1, 0))

    def test_calculate_azimuth_oceanography(self):
        self.assertEqual(0, Meteorology.calculate_azimuth_degree_oceanography(0, 1))
        self.assertEqual(90, Meteorology.calculate_azimuth_degree_oceanography(1, 0))
        self.assertEqual(180, Meteorology.calculate_azimuth_degree_oceanography(0, -1))
        self.assertEqual(270, Meteorology.calculate_azimuth_degree_oceanography(-1, 0))