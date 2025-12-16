import unittest
import logging
from Config import Config
from Chemistry import convert_column_ppb_dry_air, convert_column_ppb_with_water
from Constants import DRYAIR_MOLECULAR_MASS, GRAVITY, STANDARD_ATMOSPHERIC_PRESSURE
from numpy import average
from Source import create_source
from TROPOMI import TROPOMI_for_orbit

logger = logging.getLogger(__name__)

class Chemistry_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        Config.create_log("Chemistry_test.log")

    def __init__(self, methodName = "runTest"):
        super().__init__(methodName)

    def test_convert_column_ppb_dry_air(self):
        if not(logger is None): logger.info("Starting test")
        surface_pressure = STANDARD_ATMOSPHERIC_PRESSURE
        ppb = 1 # Concentration
        actual = convert_column_ppb_dry_air(ppb, surface_pressure)
        expected = 5.706468500102926e-06
        self.assertEqual(expected, actual)

    def test_convert_column_ppb_water_sample(self):
        '''
        Sample calculations after adding 700 H2O. 
        Arbitrary value taken from TROPOMI 09556 scan over Hail Creek
        '''
        if not(logger is None): logger.info("Starting test")
        h2o = 700 # mol m-2
        surface_pressure = STANDARD_ATMOSPHERIC_PRESSURE
        dry_atmosphere =  (surface_pressure / GRAVITY / DRYAIR_MOLECULAR_MASS * 1000.0 - h2o) # dry atmosphere mol m-2
        ppb = 1 # Concentration
        actual = convert_column_ppb_with_water(ppb, dry_atmosphere, h2o, surface_pressure)
        expected = 5.710712731413824e-06
        self.assertEqual(expected, actual)

    def test_convert_column_ppb_water(self):
        '''
        Test data taken from Hail Creek orbit 09956
        '''
        orbit = "09956"
        if not(logger is None): logger.info("Starting test")
        source = create_source("HailCreek")
        tropomi = TROPOMI_for_orbit(orbit, "010302")

        (scan, pixel, _, _) = tropomi.get_pixel_for_source(source.xy)
        ppb = tropomi.get_variable_data("methane_mixing_ratio_bias_corrected")[0,scan,pixel]
        dry_air = average(tropomi.get_variable_data("dry_air_subcolumns")[0,scan,pixel,:])
        surface_pressure = tropomi.get_variable_data("surface_pressure")[0,scan,pixel]
        h2o = tropomi.get_variable_data("water_total_column")[0, scan,pixel]
        actual = convert_column_ppb_with_water(ppb, dry_air, h2o, surface_pressure)
        expected = 0.010299725995126387 # (kg m^-2)
        self.assertEqual(expected, actual)

    def test_convert_column_ppb_water_reduces_to_dry(self):
        '''
        Verify that water reduces to dry when water content is 0 and 
        '''
        ppb = 1
        surface_pressure = STANDARD_ATMOSPHERIC_PRESSURE
        h2o = 0
        dry_air = (surface_pressure / GRAVITY / DRYAIR_MOLECULAR_MASS * 1000.0 - h2o) # dry atmosphere mol m-2
        actual = convert_column_ppb_with_water(ppb, dry_air, h2o, surface_pressure)
        expected = 5.7064685001029255e-06
        self.assertEqual(expected, actual)
