from argparse import ArgumentParser, Namespace
from ERA5 import ERA5
from matplotlib.axes import Axes
from shapely import Point

class Source:
    Sources = ["HailCreek"]
    '''
    Common properties of sources of emissions

    activity: dictionary
        Key - year of emissions. Value - Activity for this year (e.g. raw coal )
    case_name: str
        Name of the folder to keep data (e.g. HailCreek -> giving folder Data/HailCreek)
    display_name: str
        Name of source to be displayed on charts
    longitude: float
        Longitude of the source
    latitude: float
        Latitiude of the source
    '''
    def __init__(self, case_name: str, display_name: str, longitude: float, latitude: float):
        self.activity = {}
        self.case_name = case_name
        self.display_name = display_name
        self.xy = Point(longitude, latitude)
        self.gridpoint = ERA5.calculate_ERA5_NearestGridPoint(self.xy)
        self.latitude = latitude 
        self.longitude = longitude 
        self.label_latitude = latitude + 0.01
        self.label_longitude = longitude + 0.04

    def __str__(self) -> str:
        a = ["Code name: {} Display name: {}".format(self.case_name, self.display_name)]
        a.append("Longitude: {}, Latitude: {}".format(self.xy.x, self.xy.y))
        a.append("ERA5 grid point: Longitude: {}, Latitude: {}".format(self.gridpoint.x, self.gridpoint.y))
        a.append("Activity data:")
        for (k,v) in self.activity.items():
            a.append("\t{}: {}".format(k,v))
        return "\n".join(a) 

    def plot_source(self, ax: Axes):
        '''
        Default way to plot a source
        '''
        ax.plot(self.longitude, self.latitude, 'go', markersize=7)
        ax.text(self.label_longitude, self.label_latitude, self.display_name, color="black")

class Source_Mine_HailCreek(Source):    
    '''
    Hail Creek mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("HailCreek", "Hail Creek", 148.3631389, -21.490934)
    
        # "Hail Creek FactSheet 2024 - final.pdf" from Glencore website
        self.open_year = 2003
        
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity = {
            2025: 10.807512,
            2024: 9.592586,
            2023: 10.092999,
            2022: 10.783931,
            2021: 10.594483,
            2020: 9.468696,
            2019: 7.660496,
            2018: 10.203613,
            2017: 9.195752,
            2016: 10.228814,
            2015: 11.614548,
            2014: 11.911568,
            2013: 12.063330,
            2012: 12.817071,
            2011: 11.983327,
            2010: 12.554146,
        }
    
class Source_Mine_MoranbahNorth(Source):    
    '''
    MoaranbahNorth mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("MoranbahNorth", "Moranbah North", 147.948876, -21.744594)
        self.open_year = None       
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity = {
            2025: 2.269233,
            2024: 3.554822,
            2023: 5.880889,
            2022: 4.124200,
            2021: 5.849695,
            2020: 7.325862,
            2019: 7.744618,
            2018: 9.429130,
            2017: 7.785700,
            2016: 5.938807,
            2015: 7.062973,
            2014: 6.618071,
            2013: 6.169404,
            2012: 3.434443,
            2011: 5.346132,
            2010: 4.199348,
        }

def create_source(name: str) -> Source:
    '''
    Factory method for creating a source

    Parameters:
    -----------
    casename: str
        case name for source

    Returns:
    --------
    source: Source
        Source or none
    '''
    if name == "HailCreek": return Source_Mine_HailCreek()
    if name == "MoranbahNorth": return Source_Mine_MoranbahNorth()
    raise ValueError("Unknown name: {}".format(name))

def _handler_list(args: Namespace) -> None:
    source = create_source(args.source)
    print(source)

if __name__ == '__main__':
    '''
    Entry point into this module
    '''
    parser = ArgumentParser(prog="Source", description="Utility to list available sources")
    subparsers = parser.add_subparsers()

    parser_list = subparsers.add_parser("list", help="List source")
    parser_list.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_list.set_defaults(func=_handler_list)

    args = parser.parse_args()
    args.func(args)  