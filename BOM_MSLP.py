from argparse import ArgumentParser, Namespace
from Config import Config
from datetime import timedelta
from http.client import HTTPResponse
from logging import error, info, getLogger
import matplotlib.pyplot as plt
from os import path
from PIL import Image
from Source import create_source, Source
from TROPOMI import TROPOMI, TROPOMI_for_orbit
from urllib.error import URLError
from urllib.request import urlopen, Request

logger = getLogger(__name__)

def chart_MSLP(sourcename: str, orbit: str, processor: str, fileout:str) -> None:
    '''
    Chart 4 consecutive images. This is public so it can be reused by Paper.py

    Parameters
    ----------
        figure: int
            Figure counter to increment
        sourcename: str
            Identifier for source used by create_source
        orbit: str
            Orbit to identify datetime of the pass over source
        fileout:str
            Name of output or None (displays to screen)

    '''
    fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True)
    source = create_source(sourcename)
    config = Config()
    tropomi = TROPOMI_for_orbit(orbit, processor)
    [_, _, dt, _] = tropomi.get_pixel_for_source(source.xy)
    if dt is None:
        m = "Error: Orbit {} does not contain source".format(orbit)
        print(m)
        if not logger is None: logger.error(m)
        return

    dt_previous = dt - timedelta(days=1)
    year_previous = "{:04d}".format(dt_previous.year)
    month_previous  = "{:02d}".format(dt_previous.month)
    day_previous = "{:02d}".format(dt_previous.day)
    
    year_current = "{:04d}".format(dt.year)
    month_current = "{:02d}".format(dt.month)
    day_current = "{:02d}".format(dt.day)  

    file1 = path.join(Config.BoM_folder, sourcename, orbit, "IDX0102.{}{}{}1200.gif".format(year_previous, month_previous, day_previous))
    file2 = path.join(Config.BoM_folder, sourcename, orbit, "IDX0102.{}{}{}1800.gif".format(year_previous, month_previous, day_previous))
    file3 = path.join(Config.BoM_folder, sourcename, orbit, "IDX0102.{}{}{}0000.gif".format(year_current, month_current, day_current))
    file4 = path.join(Config.BoM_folder, sourcename, orbit, "IDX0102.{}{}{}0600.gif".format(year_current, month_current, day_current))

    if not path.exists(file1):
        m = "Error: Missing file: {}".format(file1)
        print(m)
        if not logger is None: logger.error(m)
        return
    if not path.exists(file2):
        m = "Error: Missing file: {}".format(file2)
        print(m)
        if not logger is None: logger.error(m)
        return
    if not path.exists(file3):
        m = "Error: Missing file: {}".format(file3)
        print(m)
        if not logger is None: logger.error(m)
        return
    if not path.exists(file4):
        m = "Error: Missing file: {}".format(file4)
        print(m)
        if not logger is None: logger.error(m)
        return
    
    img1 = Image.open(file1)
    img2 = Image.open(file2)
    img3 = Image.open(file3)
    img4 = Image.open(file4)

    axs[0,0].imshow(img1, )
    axs[0,1].imshow(img2)
    axs[1,0].imshow(img3)
    axs[1,1].imshow(img4)

    axs[0,0].axis("off")
    axs[0,1].axis("off")
    axs[1,0].axis("off")
    axs[1,1].axis("off")

    if not fileout is None:
        plt.savefig(fileout)
        print("Chart generated: {}".format(fileout))
    else:
        plt.show()

    plt.close("all")

    return

def __downloadmslp(source:str, orbit:str, year:str, month:str, day:str, hour:str):
    '''
    Download implementation
    '''
    Config.create_directory_structure_BOM(source, orbit)
    filetemplate = "IDX0102.{}{}{}{}00.gif"
    file = filetemplate.format(year, month, day, hour)
    urltemplate =  "http://www.bom.gov.au/archive/charts/{}/{}/{}"
    url = urltemplate.format(year, month, file)
    filename = path.join(Config.BoM_folder, source, orbit, file)
    if not(logger is None): info("URL: {}".format(url))

    if path.exists(filename): 
        if not(logger is None): info("File already exists file: {}".format(filename))
        return
    
    try:
        if not(logger is None): info("Downloading file: {}".format(url))
        headers = {
            "accept": "image/gif",
            "host": "www.bom.gov.au",
            "user-agent": "Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/136.0.0.0 Mobile Safari/537.36 Edg/136.0.0.0"
        }
        req = Request(url=url, headers=headers, method="GET")
        result: HTTPResponse = urlopen(req)
        bytes = result.read()
        with open(filename, "wb") as f:
            f.write(bytes)
            f.close
        print("File downloaded: {}".format(filename))
    except URLError as inst:
        if not(logger is None): error(inst.reason)
        print("File not found {}".format(filename))

def _handler_chart_MSLP(args: Namespace):
    '''
    Download 5 MSLP charts for orbit starting from 06Z on a previous day and finishing at 06Z on the image day
    Place the downloaded file in the Data,Source,Orbit,Outputs folder.
    As this is the IDX0102 product all the files will have names IDX0102.YYYYMMDDHHMM.gif
    '''
    Config.create_log("BOM_MSLP.log")
    chart_MSLP(args.source, args.orbit, args.processor, None)
    return 

def _handler_download_MSLP(args: Namespace):
    '''
    Download 5 MSLP charts for orbit starting from 06Z on a previous day and finishing at 06Z on the image day
    Place the downloaded file in the Data,Source,Orbit,Outputs folder.
    As this is the IDX0102 product all the files will have names IDX0102.YYYYMMDDHHMM.gif
    '''
    Config.create_log("BOM_MSLP.log")
    source = create_source(args.source)
    tropomi = TROPOMI_for_orbit(args.orbit, args.processor)
    [_, _, dt, _] = tropomi.get_pixel_for_source(source.xy)
    if dt is None:
        m = "Error: Orbit {} does not contain source".format(args.orbit)
        print(m)
        if not logger is None: logger.error(m)
        return
    
    dt_previous = dt - timedelta(days=1)
    year_previous = "{:04d}".format(dt_previous.year)
    month_previous  = "{:02d}".format(dt_previous.month)
    day_previous = "{:02d}".format(dt_previous.day)
    
    year_current = "{:04d}".format(dt.year)
    month_current = "{:02d}".format(dt.month)
    day_current = "{:02d}".format(dt.day)  

    __downloadmslp(args.source, args.orbit, year_previous, month_previous, day_previous, "06")
    __downloadmslp(args.source, args.orbit, year_previous, month_previous, day_previous, "12")
    __downloadmslp(args.source, args.orbit, year_previous, month_previous, day_previous, "18")
    __downloadmslp(args.source, args.orbit, year_current, month_current, day_current, "00")
    __downloadmslp(args.source, args.orbit, year_current, month_current, day_current, "06")

    return 

if __name__ == '__main__': 
    '''
    User interface for BOM MSLP download
    '''
    parser = ArgumentParser(prog="BOM_MSLP")
    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_chart = subparsers.add_parser("chart", help="Chart MSLP for orbit")
    parser_chart.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_chart.add_argument("orbit", type=str, help="The orbit number as a 5 digit. Note that the first orbit available on NCI THREDDS is 02818 and the first orbit on January 1 2020 is 11487")
    parser_chart.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart.set_defaults(func=_handler_chart_MSLP)

    parser_download = subparsers.add_parser("download", help="Download MSLP files for orbit")
    parser_download.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_download.add_argument("orbit", type=str, help="The orbit number as a 5 digit. Note that the first orbit available on NCI THREDDS is 02818 and the first orbit on January 1 2020 is 11487")
    parser_download.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_download.set_defaults(func=_handler_download_MSLP)

    args = parser.parse_args()
    
    if ("orbit" in args):
        if len(args.orbit) == 5 and args.orbit.isdigit():
            args.func(args)
            exit()
        else:
            print ("Orbit must be a five digit number or 'filtered': {}".format(args.orbit))
            exit()

    if ("processor" in args):
        if not(len(args.processor) == 6) or not(args.processor.isdigit()):
            print ("Processor must be a 6 digit number: {}".format(args.processor))
            exit()

    args.func(args)