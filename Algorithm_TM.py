from argparse import ArgumentParser, Namespace
from Algorithm_HYSPLIT import Algorithm_HYSPLIT, Algorithm_HYSPLIT_TrajectoryPoint, Algorithm_HYSPLIT_from_trajectory
from Algorithm_CSF import Algorithm_CSF, EnhancementLength, EnhancementWidth, Rotation, Transect
from AreaCalculator import AreaCalculator
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from Chemistry import convert_column_ppb_with_water
from Config import Config
from datetime import datetime, timedelta, timezone
from Geometry import Geometry
import logging
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from numpy.ma import average, count, masked, max, min
from os import listdir, path, remove
from pickle import load
from shapely import Polygon, prepare
from Source import Source, create_source
from trajectory import drawTrajectory
from TROPOMI import TROPOMI

logger = logging.getLogger(__name__)

class Algorithm_TM:

    def __init__(self, source:Source, orbit:str):
        '''
        Implementation of total mass method. It is used when the whole plume fits into the box so there is no flow through the side.
        1. The box containing plume is determined from method used by Sadavarte et all
        2. The emission time for emissions within the box is determined using HYSPLIT data
        3. Background by default is using background used by Sadavarte or can be set manually through the command line.
        '''
        self.background:float = None
        self.hours: int = None # (h)
        self.ief = None
        self.total: float = None # (kg)
        self.orbit:str = orbit
        self.q_valid_hour: float = None # (t/hr)
        self.q_valid_year: float  = None # (Gg/year)
        self.hysplit: Algorithm_HYSPLIT = None
        self.sad: Algorithm_CSF = None
        self.source = source

    def _calculate_emissions(self):
        '''
        Calculate emissions by summation of excess
        '''
        [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, _, _, scan_downwindbox]  = self.sad.tropomi.narrow_to_domain(Geometry.INTERSECTS, self.sad.downwind_box)
        pressure = self.sad.tropomi.get_variable_data("surface_pressure")[0, minscan:maxscan, minpixel:maxpixel]
        dry_air = self.sad.tropomi.get_variable_data("dry_air_subcolumns")[0,minscan:maxscan, minpixel:maxpixel,:]
        h2o = self.sad.tropomi.get_variable_data("water_total_column")[0, minscan:maxscan, minpixel:maxpixel]
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        counter = 0
        self.total = 0
        downwindbox = self.sad.downwind_box
        prepare(downwindbox)
        # For each pixel in this domain
        for s in range(scans):
            for p in range(pixels):                
                if scan_downwindbox.mask[s,p]: continue # Missing data
                v = scan_downwindbox[s,p] - self.background
                if v < 0 : continue # This is background

                pixel = Polygon([
                    (lonbox[s,p,0], latbox[s,p,0]),
                    (lonbox[s,p,1], latbox[s,p,1]),
                    (lonbox[s,p,2], latbox[s,p,2]),
                    (lonbox[s,p,3], latbox[s,p,3])
                ])
                if not(downwindbox.intersects(pixel)): continue # Not within downwind box

                intersection = Polygon(pixel.intersection(downwindbox).normalize())
                area = AreaCalculator.calculate_area(self.source, intersection)

                counter += 1
                kgm2 = convert_column_ppb_with_water(v, average(dry_air[s,p:]), h2o[s,p], pressure[s,p])
                self.total += area *  kgm2
                if not logger is None: logger.info("scan: {}, pixel: {}, area: {:.3f}, ime: {:.3f}, intersection: [{}]".format(s, p, area, area * kgm2, intersection))

        self.q_valid_hour = self.total / self.hours * 1E-3 # (t/h)
        self.q_valid_year = self.q_valid_hour * 24 * 365 * 1E-3 # (Gg/y)
        self.ief = (2 * self.q_valid_year) / (self.source.activity[2019] + self.source.activity[2020])
        return

    def chart(self, imagefile: str):
        '''
        Chart positive pixels within downwind box, and HYSPLIT 
        '''
        
        title = "TM. Orbit: {} Processor: {}\n Background {:0.3f} ppb. HYSPLIT emission puffs".format(self.sad.tropomi.orbit, self.sad.tropomi.processor_version, self.background)

        t = self.sad.tropomi_source_date
        dt_image_nearest_hour = datetime(t.year, t.month, t.day, t.hour, tzinfo=timezone.utc)
        if 30 <= t.minute: dt_image_nearest_hour += timedelta(hours=1)

        # Takes the largest box
        [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, lons, lats, scanCH4] = self.sad.tropomi.narrow_to_domain(Geometry.INTERSECTS, self.sad.downwind_box)
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        scan1 = scanCH4 - self.background
        counter = count(scan1)
        v_min = 0
        v_max = max(scan1)

        # Determine image height
        (lonmin, latmin, lonmax, latmax) = self.sad.downwind_box.bounds
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
        
        for s in range(scans):
            for p in range(pixels):
                pixelbox = Polygon([(lonbox[s,p,0], latbox[s,p,0]), (lonbox[s,p,1], latbox[s,p,1]), (lonbox[s,p,2], latbox[s,p,2]), (lonbox[s,p,3], latbox[s,p,3])])
                
                # We are dealing with multi polygon so pixels intersecting the downwind box
                # can be a point, triangle, quadrilateral or pentagon
                if self.sad.downwind_box.intersects(pixelbox) and 0 < scan1[s,p]:
                    intersection = Polygon(pixelbox.intersection(self.sad.downwind_box).normalize())
                    color = plt.cm.rainbow((scan1[s,p]-v_min)/(v_max - v_min))
                    ax.add_patch(mpatches.Polygon(intersection.exterior.coords, closed=True, facecolor=color, transform=ccrs.PlateCarree()))

        
        if (not(logger is None)): 
            logger.info("Scan count: {} min: {} max: {}".format(counter, v_min, v_max))

        # Draw value bar below chart
        sm = cm.ScalarMappable(cmap = plt.cm.rainbow, norm = mcolors.Normalize(vmin=v_min, vmax=v_max))
        cbar = fig.colorbar(sm, ax= ax, orientation="horizontal", pad=0.075)
        cbar.ax.set_xlabel("CH4 ppb")

        if not self.hysplit is None:
            points = self.hysplit.points_extent(self.sad.downwind_box, dt_image_nearest_hour)
            for p in points:
                if p.date == dt_image_nearest_hour: p.draw_point(plt, geodetic, plateCarree)

        # Draw mine marker
        ax.plot(self.sad.source.xy.x, self.sad.source.xy.y, 'go', markersize=7)
        ax.text(self.sad.source.label_longitude, self.sad.source.label_latitude, self.sad.source.display_name, color="black")

        # Draw downwind box
        drawTrajectory(plt, geodetic, plateCarree, self.sad.downwind_box.exterior.coords)

        plt.title(title)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'}

        text1 = "Background:{}(ppb), TM:{:.3f}(kg), Q: {:.3f}(t/hour) = {:.3f}(Gg/year), IEF:{:.3f}".format(self.background, self.total, self.q_valid_hour, self.q_valid_year, self.ief)
        fig.text(0.1, 0.05, text1, horizontalalignment="left", wrap=False) 
        fig.subplots_adjust(bottom=0.2)

        if not imagefile is None:
            plt.savefig(imagefile)
            print("Chart generated: {}".format(imagefile))
        else:
            plt.show()

        plt.close("all")

def Algorithm_TM_factory(sad: Algorithm_CSF, hysplit: Algorithm_HYSPLIT, background: float) -> Algorithm_TM:
        '''
        Factory method

        Parameters
        ---------
        sad: Algorithm_CSF
            Algorithm sadavarte used for downwind box determination and possibly background
        hysplit: Algorithm_HYSPLIT
            Algorithm HSYPLIT used for determination of emission time within downwind box
        background: float
            Optional background value. If None than background value from Algorithm_Sadavarte is used
        '''
        if sad is None:
            message = "Error: To run this model you need to run Algorithm Sadavarte first"
            print(message)
            if not logger is None: logger.error(message)
            return None

        if hysplit is None:
            message = "Error: To run this model you need to create Algorithm HYSPLIT first"
            print(message)
            if not logger is None: logger.error(message)
            return None

        hours = len(hysplit.points_extent(sad.downwind_box, hysplit.endtime))
        if hours is None or hours == 0:
            message = "Error: HYSPLIT trajectories do not cross downwind box. Algorithm terminates."
            print(message)
            if not logger is None: logger.error(message)
            return None

        instance = Algorithm_TM(sad.source, sad.tropomi.orbit)
        instance.hysplit = hysplit
        instance.sad = sad
        instance.background = sad.background if background is None else background
        instance.hours = hours
        instance._calculate_emissions()
        return instance

def _handler_chart(args: Namespace):
    config = Config()
    Config.create_log("Algorithm_TM_{}_{}.log".format(args.source, args.orbit))

    if args.source == "HailCreek":
        source = create_source("HailCreek")

    filename = Algorithm_CSF.get_picklename(config, source.case_name, args.orbit, args.processor)
    if not path.exists(filename):
        print("Error: Unable to opent Algorithm_CSF cache file {}. Please run the sadavart algorithm first".format(filename))
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sourcepath = path.join(Config.HYSPLIT_folder, source.case_name)
    filenames = listdir(sourcepath)
    pattern = "HYSPLIT_{}_F_".format(args.model)
    date = datetime(sad.tropomi_source_date.year, sad.tropomi_source_date.month, sad.tropomi_source_date.day, sad.tropomi_source_date.hour, tzinfo=timezone.utc)
    ah = Algorithm_HYSPLIT(source, args.model, date)
    for f in filenames:
        if not f.startswith(pattern): continue
        if not logger is None: logger.info("Processing file {}".format(f))
        ah.addTrajectories(Algorithm_HYSPLIT_from_trajectory(source, args.model, date, path.join(sourcepath, f)).trajectories)

    if 0 == len(ah.trajectories):
        if not logger is None: logger.error("No HYSPLIT trajectory files found")
        print("No HYSPLIT trajectory files found")
        return
    
    tm = Algorithm_TM_factory(sad, ah, args.background)
    tm.chart(None)

def _handler_list(args: Namespace):
    '''
    Driver for running of TM algorithm
    '''
    config = Config()
    Config.create_log("Algorithm_TM_{}_{}.log".format(args.source, args.orbit))

    if args.source == "HailCreek":
        source = create_source("HailCreek")

    filename = Algorithm_CSF.get_picklename(config, source.case_name, args.orbit, args.processor)
    if not path.exists(filename):
        print("Error: Unable to opent Algorithm_CSF cache file {}. Please run the sadavart algorithm first".format(filename))
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sourcepath = path.join(Config.HYSPLIT_folder, source.case_name)
    filenames = listdir(sourcepath)
    pattern = "HYSPLIT_{}_F_".format(args.model)
    date = datetime(sad.tropomi_source_date.year, sad.tropomi_source_date.month, sad.tropomi_source_date.day, sad.tropomi_source_date.hour, tzinfo=timezone.utc)
    ah = Algorithm_HYSPLIT(source, args.model, date)
    for f in filenames:
        if not f.startswith(pattern): continue
        if not logger is None: logger.info("Processing file {}".format(f))
        ah.addTrajectories(Algorithm_HYSPLIT_from_trajectory(source, args.model, date, path.join(sourcepath, f)).trajectories)

    if 0 == len(ah.trajectories):
        if not logger is None: logger.error("No HYSPLIT trajectory files found")
        print("No HYSPLIT trajectory files found")
        return
    
    tm = Algorithm_TM_factory(sad, ah, args.background)    
    print("background:{}(ppb), tm:{:.3f}(kg), Q: {:.3f}(t/hour) = {:.3f}(Gg/year), IEF:{:.3f}".format(tm.background, tm.total, tm.q_valid_hour, tm.q_valid_year, tm.ief))

if __name__ == '__main__':
    parser = ArgumentParser(prog="Algorithm_TM")
    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_chart = subparsers.add_parser("chart", help="Chart model")
    parser_chart.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions", )
    parser_chart.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 09956")
    parser_chart.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart.add_argument("model", type=str, choices=["CDC", "GDAS", "GFSQ"])
    parser_chart.add_argument("background", type=float, help="Background value")
    parser_chart.set_defaults(func=_handler_chart)

    parser_list = subparsers.add_parser("list", help="List model results")
    parser_list.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions", )
    parser_list.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 09956")
    parser_list.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_list.add_argument("model", type=str, choices=["CDC", "GDAS", "GFSQ"])
    parser_list.add_argument("background", type=float, help="Background value")
    parser_list.set_defaults(func=_handler_list)

    args = parser.parse_args()
    args.func(args)
