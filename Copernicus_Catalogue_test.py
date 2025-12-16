from Config import Config
from Copernicus_Catalogue import Copernicus_Catalogue
import logging
import time
import unittest

logger = logging.getLogger(__name__)

class Copernicus_Catalogue_test(unittest.TestCase):
    def __init__(self, methodName = "runTest"):
        super().__init__(methodName)

    @classmethod
    def setUpClass(cls):
        Config.create_log("Copernicus_Catalogue_test.log")

    @classmethod
    def tearDownClass(cls):
        pass

    def test_authentication(self):
        '''
        Test for orbit Hail Creek 09956
        '''
        config = Config()
        username = input("Username: ")
        password = input("Password: ")
        catalogue = Copernicus_Catalogue(username, password, config)
        r = catalogue._authenticate()
        self.assertTrue(r)
        r1 = catalogue._refreshtoken()
        self.assertTrue(r1)
        time.sleep(11*60)
        r2 = catalogue._authenticate()
        self.assertTrue(r2)

