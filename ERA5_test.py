from Config import Config
from ERA5 import ERA5, ERA5_PressureLevels, ERA5_SingleLevel
import logging
from os import path
from shapely import Point
from Source import create_source
import unittest

logger = logging.getLogger(__name__)

class ERA5_test(unittest.TestCase):

    def __init__(self, methodName = "runTest"):
        super().__init__(methodName)

    @classmethod
    def setUpClass(cls):
        Config.create_log("ERA5_Test.log")
    
    @classmethod
    def tearDownClass(cls):
        pass

    def test_calculate_ERA5_NearestGridPoint(self):
        source = create_source("HailCreek")
        actual = ERA5.calculate_ERA5_NearestGridPoint(source.xy)
        expected = Point(148.25, -21.5)
        self.assertEqual(actual.x, expected.x)
        self.assertEqual(actual.y, expected.y)

    def test_download_SingleLevel(self):
        source = create_source("HailCreek")
        orbit = "09956"
        Config.create_directory_structure_ERA5(source.case_name)

        target = ERA5_SingleLevel.get_filepath("CSF", source.case_name, orbit)
        if path.exists(target):
            self.assertTrue(True)
        else:
            request = {
                "product_type": ["reanalysis"],
                "variable": ["boundary_layer_height", "geopotential", "surface_pressure", "10m_u_component_of_wind", "10m_v_component_of_wind"],
                "year": ["2019"],
                "month": ["09"],
                "day": ["15"],
                "time": ["04:00"],
                "data_format": "netcdf",
                "download_format": "unarchived",
                "area": [-21.5, 148.25, -21.5, 148.25]
            }
            e5 = ERA5_SingleLevel.download(request, target)
            self.assertTrue(path.exists(target))

    def test_download_PressureLevels(self):
        
        # start Test getting variable. Read u at pressure level 1000 hPa
        source = create_source("HailCreek")
        orbit = "09956"
        Config.create_directory_structure_ERA5(source.case_name)
       
        target = ERA5_PressureLevels.get_filepath("CSF", source.case_name, orbit)
        if path.exists(target):
            self.assertTrue(True)
        else:
            request = {
                "product_type": ["reanalysis"],
                "variable": [
                    "geopotential",
                    "u_component_of_wind",
                    "v_component_of_wind"
                ],
                "year": ["2019"],
                "month": ["09"],
                "day": ["15"],
                "time": ["04:00"],
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
                "area": [-21.5, 148.25, -21.5, 148.25]
            }
            e5 = ERA5_PressureLevels.download(request, target)
            self.assertTrue(path.exists(target))
