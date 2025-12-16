from argparse import ArgumentParser, Namespace
from Algorithm_CSF import Algorithm_CSF, EnhancementLength, EnhancementWidth, Rotation, Transect
from AreaCalculator import AreaCalculator
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from Chemistry import convert_column_ppb_with_water
from Config import Config
from datetime import datetime, timedelta, timezone
from Geometry import Geometry
import logging
from math import cos, pi, sqrt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from numpy.ma import average, count, max
from os import path, remove
from pickle import load
from shapely import Polygon, prepare
from Source import Source, create_source
from trajectory import drawTrajectory

logger = logging.getLogger(__name__)

class Algorithm_IME:

    def __init__(self, sad:Algorithm_CSF, transect:int):
        '''
        Implementation of integrated mass enhancement method
        See Varon. et all "Quantifying methane point sources from fine-scale satellite observations of atmospheric methane plumes"
        Atmos. Meas. Tech., 11, 5673 - 5686, 2018 https://doi.org/10.5194/amt-11-5673-2018

        1. The box containing plume is determined from method used by Sadavarte et all
        2. The emission time for emissions within the box is determined distance from transect to the source divided by wind speed
        3. Background by default is using background used by Sadavarte or can be set manually through the command line.
        '''
        self.background:float = None
        self.downwindbox: Polygon = None
        self.ief: float = None
        self.ime: float = None # (kg)
        self.length: float = None # (m)
        self.orbit:str = sad.tropomi.orbit
        self.q_valid_hour: float = None # (t/hr)
        self.q_valid_year: float  = None # (Gg/year)
        self.residence_hours: float = None # (h)
        self.sad: Algorithm_CSF = sad
        self.source = sad.source
        direction_correction = cos(pi * (self.sad.downwind_box_azimuth_adjusted - self.sad.downwind_box_azimuth_initial) / 180.0)
        self.windspeed = sqrt(self.sad.u * self.sad.u + self.sad.v * self.sad.v) * direction_correction
        self.transectid = transect
        t = self.sad.transects[self.transectid - 4]
        if not t.isvalid(): 
            m = "Error: Transect {} is not valid.".format(self.transectid)
            if not logger is None: logger.error(m)
            print(m)
            return

        # Calculate 0 line of the box
        psource = [self.source.xy.x, self.source.xy.y]
        [[p1_lon, p1_lat, _]] = Algorithm_CSF.g.direct(points=psource, azimuths=self.sad.downwind_box_azimuth_adjusted + 90, distances=self.sad.downwind_box_half_width)
        [[p2_lon, p2_lat, _]] = Algorithm_CSF.g.direct(points=psource, azimuths=self.sad.downwind_box_azimuth_adjusted - 90, distances=self.sad.downwind_box_half_width)
        p1 = [p1_lon, p1_lat]
        p2 = [p2_lon, p2_lat]
        
        # Calculate mid point of transect
        p3 = [t.line.coords.xy[0][0], t.line.coords.xy[1][0]]
        p4 = [t.line.coords.xy[0][1], t.line.coords.xy[1][1]]
        [[pm_lon, pm_lat, _]] = Algorithm_CSF.g.direct(points=p3, azimuths=self.sad.downwind_box_azimuth_adjusted + 90, distances=self.sad.downwind_box_half_width)

        # Calculate distance to the source
        dist_inverse = Algorithm_CSF.g.inverse([(pm_lon, pm_lat)], [(self.source.longitude, self.source.latitude)])
        self.length = dist_inverse[0][0]
        # Calculate residence time in hours
        self.residence_hours = self.length / self.windspeed / 3600        
        # Sort points into a squence
        pur = None
        pul = None
        pbl = None
        pbr = None
        pall = [p1, p2, p3, p4]

        # Sort by longitude ascending
        pall.sort(key=lambda item: item[0])
        if pall[0][1] < pall[1][1] : 
            pbl = pall[0]
            pul = pall[1]
        else:
            pbl = pall[1]
            pul = pall[0]

        if pall[2][1] < pall[3][1] : 
            pbr = pall[2]
            pur = pall[3]
        else:
            pbr = pall[3]
            pur = pall[2]

        self.downwindbox = Polygon([(pul[0], pul[1]), (pbl[0], pbl[1]), (pbr[0], pbr[1]), (pur[0], pur[1])])
        prepare(self.downwindbox)

    def _calculate_emissions(self):
        '''
        Calculate emissions by summation of excess. Time is determined as ratio of distance between tranect and mine and wind speed and .
        '''
        if self.downwindbox is None: return
        if not logger is None: logger.info("IME. Orbit: {} Processor {} Background {:0.3f} ppb.".format(self.sad.tropomi.orbit, self.sad.tropomi.processor_version, self.background))

        # Calculate integrated mass within a box
        [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, _, _, scan_downwindbox]  = self.sad.tropomi.narrow_to_domain(Geometry.INTERSECTS, self.downwindbox)
        pressure = self.sad.tropomi.get_variable_data("surface_pressure")[0, minscan:maxscan, minpixel:maxpixel]
        dry_air = self.sad.tropomi.get_variable_data("dry_air_subcolumns")[0,minscan:maxscan, minpixel:maxpixel,:]
        h2o = self.sad.tropomi.get_variable_data("water_total_column")[0, minscan:maxscan, minpixel:maxpixel]
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        counter = 0
        self.ime = 0
        # For each pixel in this domain
        for s in range(scans):
            for p in range(pixels):                
                if scan_downwindbox.mask[s,p]: continue # Missing data
                v = scan_downwindbox[s,p] - self.background
                if v < 0 : continue # This is background

                bl = [lonbox[s,p,0], latbox[s,p,0]]
                ul = [lonbox[s,p,1], latbox[s,p,1]]
                ur = [lonbox[s,p,2], latbox[s,p,2]]
                br = [lonbox[s,p,3], latbox[s,p,3]]

                pixel = Polygon([ul,ur, br, bl])
                if not(self.downwindbox.intersects(pixel)): continue 

                intersection = Polygon(pixel.intersection(self.downwindbox).normalize())
                area = AreaCalculator.calculate_area(self.source, intersection)

                kgm2 = convert_column_ppb_with_water(v, average(dry_air[s,p:]), h2o[s,p], pressure[s,p])
                counter += 1
                self.ime += area *  kgm2
                if not logger is None: logger.info("scan: {}, pixel: {}, area: {:.3f}, ime: {:.3f}, intersection: [{}]".format(s, p, area, area * kgm2, intersection))

        self.q_valid_hour = self.ime / self.residence_hours * 1E-3 # (t/h)
        self.q_valid_year = self.q_valid_hour * 24 * 365 * 1E-3 # (Gg/y)
        self.ief = (2 * self.q_valid_year) / (self.source.activity[2019] + self.source.activity[2020])
        return

    def chart(self, imagefile: str):
        '''
        Chart positive pixels within downwind box 
        '''
        if self.downwindbox is None: return        
        title = "IME. Orbit: {} Processor {} \n Background {:0.3f} ppb.".format(self.sad.tropomi.orbit, self.sad.tropomi.processor_version, self.background)

        t = self.sad.tropomi_source_date
        dt_image_nearest_hour = datetime(t.year, t.month, t.day, t.hour, tzinfo=timezone.utc)
        if 30 <= t.minute: dt_image_nearest_hour += timedelta(hours=1)

        # Takes the largest box
        [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, _, _, scanCH4] = self.sad.tropomi.narrow_to_domain(Geometry.INTERSECTS, self.downwindbox)
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        scan1 = scanCH4 - self.background

        # Determine image height
        (lonmin, latmin, lonmax, latmax) = self.downwindbox.bounds
        lonwidth = lonmax - lonmin
        latheight = latmax - latmin
        mll = lonwidth
        if mll < latheight: mll = latheight
        fracy = int(10 * latheight/mll)
        if fracy < 5: fracy = 5

        plateCarree = ccrs.PlateCarree()
        geodetic = ccrs.Geodetic()
        fig = plt.figure(figsize=(10,fracy))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extents=[lonmin-0.1, lonmax+0.1, latmin-0.1, latmax+0.1], crs=plateCarree)
        ax.coastlines() 

        v_min = 0
        v_max = max(scan1)
        counter = count(scan1)       

        for s in range(scans):
            for p in range(pixels):
                pixelbox = Polygon([(lonbox[s,p,0], latbox[s,p,0]), (lonbox[s,p,1], latbox[s,p,1]), (lonbox[s,p,2], latbox[s,p,2]), (lonbox[s,p,3], latbox[s,p,3])])
                
                # We are dealing with multi polygon so pixels intersecting the downwind box
                # can be a point, triangle, quadrilateral or pentagon
                if self.downwindbox.intersects(pixelbox) and 0 < scan1[s,p]:
                    intersection = Polygon(pixelbox.intersection(self.downwindbox).normalize())
                    color = plt.cm.rainbow((scan1[s,p]-v_min)/(v_max - v_min))
                    ax.add_patch(mpatches.Polygon(intersection.exterior.coords, closed=True, facecolor=color, transform=ccrs.PlateCarree()))

        if (not(logger is None)): 
            logger.info("Scan count: {} min: {} max: {}".format(counter, v_min, v_max))
       
        # Draw value bar below chart
        sm = cm.ScalarMappable(cmap = plt.cm.rainbow, norm = mcolors.Normalize(vmin=v_min, vmax=v_max))
        cbar = fig.colorbar(sm, ax= ax, orientation="horizontal", pad=0.075)
        cbar.ax.set_xlabel("CH4 ppb")

        # Draw mine marker        
        ax.plot(self.sad.source.xy.x, self.sad.source.xy.y, 'go', markersize=7)
        ax.text(self.sad.source.label_longitude, self.sad.source.label_latitude, self.sad.source.display_name, color="black")

        # Draw downwind box
        drawTrajectory(plt, geodetic, plateCarree, self.downwindbox.exterior.coords)

        plt.title(title)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'}
        
        text1 = "Transect ID:{}, Distance to source: {:.3f}(km) ".format(self.transectid, self.length * 0.001)
        text2 = "Wind:{:.3f} (km/h), Residence time : {:.3f}(h)".format(self.windspeed * 3.6, self.residence_hours)
        text3 = "Background:{}(ppb), IME:{:.3f}(kg), Q: {:.3f}(t/hour) = {:.3f}(Gg/year), IEF:{:.3f}".format(self.background, self.ime, self.q_valid_hour, self.q_valid_year, self.ief)

        fig.text(0.1, 0.1, text1, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.075, text2, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.05, text3, horizontalalignment="left", wrap=False) 
        fig.subplots_adjust(bottom=0.2)

        if not imagefile is None:
            plt.savefig(imagefile)
            print("Chart generated: {}".format(imagefile))
        else:
            plt.show()

        plt.close("all")


def Algorithm_IME_factory(sad: Algorithm_CSF, background: float, transectid: int) -> Algorithm_IME:
        '''
        Factory method

        Parameters
        ---------
        sad: Algorithm_CSF
            Algorithm sadavarte used for downwind box determination and possibly background
        background: float
            Optional background value. If None than background value from Algorithm_Sadavarte is used
        '''
        if sad is None:
            message = "Error: To run this model you need to run Algorithm Sadavarte first"
            print(message)
            if not logger is None: logger.error(message)
            return None

        instance = Algorithm_IME(sad,  transectid)
        instance.sad = sad
        instance.background = background
        instance._calculate_emissions()
        return instance

def _handler_chart(args: Namespace):
    config = Config()
    Config.create_log("Algorithm_IME_{}_{}.log".format(args.source, args.orbit))

    if args.transectid < 4:
        print("Error: Transect id must be at least 4 {}.".format(args.transectid))
        return
    if 15 < args.transectid:
        print("Error: Transect id must be at most 15 {}.".format(args.transectid))
        return

    if args.source == "HailCreek":
        source = create_source("HailCreek")

    filename = Algorithm_CSF.get_picklename(config, source.case_name, args.orbit, args.processor)
    if not path.exists(filename):
        print("Error: Unable to opent Algorithm_CSF cache file {}. Please run the sadavart algorithm first".format(filename))
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()
    
    ime = Algorithm_IME_factory(sad, args.background if args.background else sad.background, args.transectid)
    ime.chart(None)

def _handler_list(args: Namespace):
    '''
    Driver for running of IME algorithm
    '''
    config = Config()
    Config.create_log("Algorithm_IME_{}_{}.log".format(args.source, args.orbit))

    if args.transectid < 4:
        print("Error: Transect id must be at least 4 {}.".format(args.transectid))
        return
    if 15 < args.transectid:
        print("Error: Transect id must be at most 15 {}.".format(args.transectid))
        return
    
    if args.source == "HailCreek":
        source = create_source("HailCreek")

    filename = Algorithm_CSF.get_picklename(config, source.case_name, args.orbit, args.processor)
    if not path.exists(filename):
        print("Error: Unable to opent Algorithm_CSF cache file {}. Please run the sadavart algorithm first".format(filename))
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()
    
    ime = Algorithm_IME_factory(sad, args.background if args.background else sad.background, args.transectid)  
    if ime.length is None: return
    print("Transect ID:{}, Distance to source: {:.3f}(km) ".format(ime.transectid, ime.length * 0.001))    
    print("Wind:{:.3f} (km/h), Residence time : {:.3f}(h)".format(ime.windspeed * 3.6, ime.residence_hours))    
    print("Background:{}(ppb), IME:{:.3f}(kg), Q: {:.3f}(t/hour) = {:.3f}(Gg/year), IEF:{:.3f}".format(ime.background, ime.ime, ime.q_valid_hour, ime.q_valid_year, ime.ief))

if __name__ == '__main__':
    parser = ArgumentParser(prog="Algorithm_IME")
    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_chart = subparsers.add_parser("chart", help="Chart model")
    parser_chart.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions", )
    parser_chart.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 09956")
    parser_chart.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart.add_argument("transectid", type=int, help="Transect number between 4 and 15")
    parser_chart.add_argument("-b", "--background", type=float, help="Background value. If not specified use CSF background")
    parser_chart.set_defaults(func=_handler_chart)

    parser_list = subparsers.add_parser("list", help="List model results")
    parser_list.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions", )
    parser_list.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 09956")
    parser_list.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_list.add_argument("transectid", type=int, help="Transect number between 4 and 15")
    parser_list.add_argument("-b", "--background", type=float, help="Background value. If not specified use CSF background")
    parser_list.set_defaults(func=_handler_list)
            
    args = parser.parse_args()
    if ("orbit" in args):
        if not(len(args.orbit) == 5) or not(args.orbit.isdigit()):
            print ("Orbit must be a five digit number: {}".format(args.orbit))
            exit()

    if ("processor" in args):
        if not(len(args.processor) == 6) or not(args.processor.isdigit()):
            print ("Processor must be a 6 digit number: {}".format(args.processor))
            exit()

    args.func(args)
