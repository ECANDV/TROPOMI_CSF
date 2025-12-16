from argparse import ArgumentParser
import boto3
from Config import Config
import logging
from os import path
from pickle import load
from TROPOMI import TROPOMI
from TROPOMI_Filter import TROPOMI_Filter

logger = logging.getLogger(__name__)

class Copernicus_S3:
    '''
    Download TROPOMI files from Copernicus EU S3
    '''
    def __init__(self, config: Config, aws_access_key_id: str, aws_secret_access_key: str):
        self.config = config
        self.session = boto3.session.Session()
        self.s3 = boto3.resource("s3",
            endpoint_url='https://eodata.dataspace.copernicus.eu',
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
            region_name='default'
        )
        pass
        
    def _download_single(self, fin: str) -> None:
        '''
        Download S3 url for files containing domain specified in Config file
        '''
        bucket = self.s3.Bucket("eodata")
        f = fin[-86:]
        processor = TROPOMI.get_processor_version(f)
        Config.create_directory_structure_TROPOMI(processor)
        fout = path.join(Config.TROPOMI_folder, processor, f)
        if path.exists(fout):
            if not logger is None: logger.info(f"Target file alread exists: {fout}")
            return

        try:
            bucket.download_file(fin, fout)
            if not logger is None: logger.info(f"Success: {fout}")
        except:
            if not logger is None: logger.error(f"Error: {fout}")

    def _download_filter(self) -> None:
        '''
        Download S3 url for files containing domain specified in Config file
        '''
        if not path.exists(TROPOMI_Filter.pickle):
            m = "Error: Please execute python '.\\TROPOMI_Filter run contains' command first"
            if not logger is None: logger.error(m)
            print(m)
            return
        
        with open(TROPOMI_Filter.pickle, "rb") as f:
            filter: TROPOMI_Filter = load(f)
            f.close()

        for k in filter.success_contains_s3_key:
            self._download_single(k)

def _handler_download_filter(args):
    config = Config()
    Config.create_log("Copernicus_s3.log")
    cs3 = Copernicus_S3(config, args.aws_access_key_id, args.aws_secret_access_key)
    cs3._download_filter()
    

def _handler_download_single(args):
    config = Config()
    Config.create_log("Copernicus_s3.log")
    cs3 = Copernicus_S3(config, args.aws_access_key_id, args.aws_secret_access_key)
    cs3._download_single(args.aws_key)

if __name__ == '__main__':
    parser = ArgumentParser(prog="Copernicus_S3")
    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_single = subparsers.add_parser("single", help="Download a single file")
    parser_single.add_argument("aws_access_key_id", type=str, help="Copernicus Dataspace aws_access_key_id")
    parser_single.add_argument("aws_secret_access_key", type=str, help="Copernicus Dataspace aws_secret_access_key")
    parser_single.add_argument("aws_key", type=str, help="AWS s3 key in format: Sentinel-5P/TROPOMI/L2__CH4___/2019/09/15/S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc")
    parser_single.set_defaults(func=_handler_download_single)

    parser_single = subparsers.add_parser("filter", help="Download all files selected by contain filter")
    parser_single.add_argument("aws_access_key_id", type=str, help="Copernicus Dataspace aws_access_key_id")
    parser_single.add_argument("aws_secret_access_key", type=str, help="Copernicus Dataspace aws_secret_access_key")
    parser_single.set_defaults(func=_handler_download_filter)

    args = parser.parse_args()
    args.func(args)
    
