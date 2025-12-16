from argparse import ArgumentParser
from datetime import datetime, timedelta, timezone
from Config import Config
import json
import logging
import requests

logger = logging.getLogger(__name__)

class Copernicus_Catalogue:
    url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
    headers = { "Content-Type": "application/x-www-form-urlencoded" }
    '''
    Search for TROPOMI files in Copernicus EU catalogue
    '''
    def __init__(self, username: str, password: str, config: Config):
        self.username = username
        self.password = password
        self.config = config
        self.access_token = None
        self.access_time = None
        self.expires_in = None
        self.refresh_expires_in = None
        self.refresh_token = None
        
        pass
        
    def _authenticate(self) -> bool:
        '''
        Use pasword/username authentication to obtain access token
        '''
        
        n = datetime.now(tz=timezone.utc)

        if not self.access_time is None and not self.expires_in is None:
            e_auth = self.access_time + timedelta(seconds=self.expires_in - 10)
            # Authentication is still valid for at least anohter 10 seconds. Return true
            if n < e_auth: return True
        
            # Authentication is not valid but refresh token is valid
            elif not self.refresh_token is None and not self.refresh_expires_in is None: 
                e_refresh = self.access_time + timedelta(seconds=self.refresh_expires_in - 10)
                if n < e_refresh and self._refreshtoken(): return True

        # Reset 
        self.access_token = None
        self.access_time = None
        self.expires_in = None
        self.refresh_expires_in = None
        self.refresh_token = None        
        
        data={
            "username": self.username,
            "password": self.password,
            "grant_type" : "password",
            "client_id": "cdse-public"
        }
        response = requests.post(Copernicus_Catalogue.url, headers=Copernicus_Catalogue.headers, data=data)
        if response is None:
            print("Error: Response is None")
            if not logger is None: logger.error("Error: Response is None")
            return False
        if not response.status_code == 200:
            print("Error: Status code: {}".format(response.status_code))
            print(response.text)
            if not logger is None: logger.error("Error: Status code: {}".format(response.status_code))
            if not logger is None: logger.error(response.text)
            return False
        j = json.loads(response.text)
        self.access_token = j["access_token"] 
        self.expires_in = float(j["expires_in"])
        self.access_time = datetime.now(tz=timezone.utc)
        self.refresh_expires_in = float(j["refresh_expires_in"])
        self.refresh_token = j["refresh_token"]        
        if not logger is None:
            logger.info(self.access_time)
            logger.info(self.access_token)
            logger.info(self.expires_in)
            logger.info(self.refresh_expires_in)
            logger.info(self.refresh_token)

        return True

    def _refreshtoken(self) -> bool:
        
        data={
            "grant_type" : "refresh_token",
            "refresh_token": self.refresh_token,
            "client_id": "cdse-public"
        }
        # Reset 
        self.access_token = None
        self.access_time = None
        self.expires_in = None
        self.refresh_expires_in = None
        self.refresh_token = None        

        response = requests.post(Copernicus_Catalogue.url, headers=Copernicus_Catalogue.headers, data=data)
        if response is None:
            print("Error: Response is None")
            if not logger is None: logger.error("Error: Response is None")
            return False
        if not response.status_code == 200:
            print("Error: Status code: {}".format(response.status_code))
            print(response.text)
            if not logger is None: logger.error("Error: Status code: {}".format(response.status_code))
            if not logger is None: logger.error(response.text)
            return False

        j = json.loads(response.text)
        self.access_time = datetime.now(tz=timezone.utc)
        self.access_token = j["access_token"] 
        self.expires_in = float(j["expires_in"])
        self.refresh_expires_in = float(j["refresh_expires_in"])
        self.refresh_token = j["refresh_token"]
        if not logger is None:
            logger.info("access_time: {}".format(self.access_time))
            logger.info("access_token: {}".format(self.access_token))
            logger.info("expires_in: {}".format(self.expires_in))
            logger.info("refresh_expires_in: {}".format(self.refresh_expires_in))
            logger.info("refresh_token: {}".format(self.refresh_token))
        return True

    def get_s3_keys(self, date: str) -> list[str]:
        '''
        List S3 keys for files containing domain specified in Config file for date

        Parameters
        ----------
        date: str
            Date in YYYYMMDD format
        '''
        results = []
        if not self._authenticate():
            m = "Error: Failure to authenticate. Please check you password username"
            print(m)
            if not logger is None: logger.error(m)
            # Want to stop iterations
            raise ConnectionRefusedError() 
        year = date[0:4]
        month = date[4:6]
        day = date[6:8]
        url = "https://sh.dataspace.copernicus.eu/api/v1/catalog/1.0.0/search" 
        headers = {   
            "Content-Type": "application/json",   
            "Authorization": "Bearer " + self.access_token 
        } 
        data = {   
            "collections": ["sentinel-5p-l2"],   
            "datetime": "{y}-{m}-{d}T00:00:00Z/{y}-{m}-{d}T23:59:59Z".format(y=year, m=month, d=day),   
            "bbox": [self.config.Algorithm_CSF_domain_xmin,self.config.Algorithm_CSF_domain_ymin,self.config.Algorithm_CSF_domain_xmax,self.config.Algorithm_CSF_domain_ymax],   
            "limit": 10,   
            "filter": {
                "op": "=",
                "args": [{"property": "s5p:type"},"CH4"]
            }, 
            "filter-lang": "cql2-json" }  

        response = requests.post(url, headers=headers, json=data)
        if response is None:
            print("Error: Response is None")
            return results
        if not response.status_code == 200:
            print("Error: Status code: {}".format(response.status_code))
            print(response.text)
            return results   
        
        s = json.loads(response.text)
        features = s["features"]
        for f in features:
            if "id" in f.keys(): 
                # Strip service and bucket s3://EODATA/ from response
                s = "{}/{}".format(f["assets"]["data"]["href"], f["id"])[12:]
                results.append(s)
                print(s)
        return results

def _handler_list(args):
    config = Config()
    cc = Copernicus_Catalogue(args.username, args.password, config)
    cc.get_s3_keys(args.date)

if __name__ == '__main__':
    parser = ArgumentParser(prog="Copernicus_Catalogue")
    parser.add_argument("date", type=str, help="Date in format YYYYMMDD")
    parser.add_argument("username", type=str, help="Copernicus EU username")
    parser.add_argument("password", type=str, help="Copernicus EU password")
    parser.set_defaults(func=_handler_list)
    args = parser.parse_args()
    args.func(args)
    
