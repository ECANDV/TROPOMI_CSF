from argparse import ArgumentParser, Namespace
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from Config import Config
from datetime import datetime, timedelta, timezone
from Geometry import Geometry
import logging
from matplotlib.colors import BoundaryNorm, ListedColormap
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from numpy import int8, zeros
from numpy.ma import average, max, min
from os import listdir, path, remove
from shapely import box, Point, Polygon, prepare
from Source import Source, create_source
from TROPOMI import TROPOMI, TROPOMI_for_orbit

logger = logging.getLogger(__name__)

class Algorithm_HYSPLIT_TrajectoryPoint:
    '''
    A point in HYSPLIT trajectory. Taken from 
    I6 - trajectory number
    I6 - meteorological grid number
    5I6 - time of point: year month day hour minute
    I6 - forecast hour at point
    F8.1 - age of the trajectory in hours
    2F9.3 - position latitude and longitude
    F8.1 - position height in meters above ground
    n1XF8.1 - n diagnostic output variables (1st output is always pressure)
    '''
    def __init__(self, 
        gridnumber:int, 
        date: datetime, 
        forecasthour: int, 
        ageoftrajectory: float, 
        latitude: float, 
        longitude: float, 
        height: float, 
        diagnostics:list[float]
        ):
        self.gridnumber = gridnumber
        self.date = date
        self.forecasthour = forecasthour
        self.ageoftrajectory = ageoftrajectory
        self.latitude = latitude
        self.longitude = longitude
        self.height = height
        self.diagnostics = diagnostics
        self.pressure = self.diagnostics[0]

    def __eq__(self, value) -> bool:
        '''
        Equality of trajectory points
        '''
        if not isinstance(value, Algorithm_HYSPLIT_TrajectoryPoint): return False
        if self.gridnumber != value.gridnumber: return False
        if self.date != value.date: return False
        if self.forecasthour != value.forecasthour: return False
        if self.ageoftrajectory != value.ageoftrajectory: return False
        if self.latitude != value.latitude: return False
        if self.longitude != value.longitude: return False
        if self.height != value.height: return False
        if len(self.diagnostics) != len(value.diagnostics): return False
        for i in range(len(self.diagnostics)):
            if self.diagnostics[i] != value.diagnostics[i]: return False
        return True

    def __str__(self):
        msgs = []
        msgs.append("Point: Date: {}, Longitude: {}, Latitude: {}".format(self.date, self.longitude, self.latitude))
        msgs.append("\tGrid: {}".format(self.gridnumber))
        msgs.append("\tForecast hour: {}".format(self.forecasthour))
        msgs.append("\tPuff age: {}".format(self.ageoftrajectory))
        msgs.append("\tHeight: {}".format(self.height))
        msgs.append("\tPressure: {}".format(self.pressure))
        return "\n".join(msgs)
        
    def draw_point(self, plt, geodetic, plateCarree) -> None:
        '''
        Draw point if it belongs to the image
        '''
        s_lon, s_lat = plateCarree.transform_point(self.longitude, self.latitude, geodetic)
        m = 'gX' if 0 <= self.ageoftrajectory else 'gP'
        plt.plot(s_lon, s_lat, m, markersize=10, transform=plateCarree)

class Algorithm_HYSPLIT_Trajectory:
    '''
    HYSPLIT trajectory
    '''
    def __init__(self, fixed_latitude: float, fixed_longitude: float, fixed_height: float):
        '''
        Definition of trajectory as used by HYSPLIT

        Parameters
        ----------
        fixed_latitude - Fixed point of trajectory latitiude
        fixed_longitude - Fixed point of trajectory longitude
        fixed_height - Fixed point height (m) above the ground

        '''
        self.fixed_latitude = fixed_latitude
        self.fixed_longitude = fixed_longitude
        self.fixed_height = fixed_height
        self.points: dict[datetime, Algorithm_HYSPLIT_TrajectoryPoint] = dict() # key: datetime value Trajectory points

    def __eq__(self, value) -> bool:
        '''
        Trajectory equality
        '''
        if not isinstance(value, Algorithm_HYSPLIT_Trajectory): return False
        v1: Algorithm_HYSPLIT_Trajectory = value
        if not self.fixed_latitude == v1.fixed_latitude : return False
        if not self.fixed_longitude == v1.fixed_longitude : return False
        if not self.fixed_height == v1.fixed_height : return False

        lp = len(self.points)
        vp = len(v1.points)
        if lp == 0 and vp == 0: return True
        if not lp == vp: return False
        for i in self.points.keys():
            if not i in v1.points: return False
            if not self.points[i] == v1.points[i]: return False
        return True

    def __str__(self) -> str:
        msgs = []
        msgs.append("Trajectory: Points")
        for p in self.points:
            msgs.append(p.__str__())
        return "\n".join(msgs)

    def draw_trajectory(self, plt, geodetic, plateCarree, end_datetime: datetime) -> None:
        p_p = None
        p_n = None
        for p in self.points.keys():
            p_n = self.points[p]
            if not p_p is None:
                s_lon, s_lat = plateCarree.transform_point(p_p.longitude, p_p.latitude, geodetic)
                e_lon, e_lat = plateCarree.transform_point(p_n.longitude, p_n.latitude, geodetic)

                m ="X" if p_p.ageoftrajectory < 0 else "P"
                plt.plot([s_lon, e_lon], [s_lat, e_lat], color="black", markersize=10,  markerfacecolor="green", linewidth=1, marker=m, transform=plateCarree)
            
            p_p = p_n

    def end_datetime(self) -> datetime:
        dt = None
        for k in self.points.keys():
            if dt is None or dt < k: dt = k
        return dt

    def end_latitude(self) -> float:
        dt = self.end_datetime()
        if dt is None: return None
        return self.points[dt].latitude

    def end_longitude(self) -> float:
        dt = self.end_datetime()
        if dt is None: return None
        return self.points[dt].longitude

    def start_datetime(self) -> datetime:
        dt = None
        for k in self.points.keys():
            if dt is None or k < dt: dt = k
        return dt

    def start_latitude(self) -> float:
        dt = self.start_datetime()
        if dt is None: return None
        return self.points[dt].latitude

    def start_longitude(self) -> float:
        dt = self.start_datetime()
        if dt is None: return None
        return self.points[dt].longitude

class Algorithm_HYSPLIT:
    '''
    HYSPLIT model containing a list of trajectories from HYSPLIT.
    '''
        
    def __init__(self, source:Source, modelstr:str, date_end:datetime):
        '''
        Representation of HYSPLIT trajectory modelling output
        '''
        self.background = None
        self.diagnostics = []
        self.direction = None
        self.endtime = date_end # datetime(int(datestr[0:4]), int(datestr[4:6]), int(datestr[6:8]), int(datestr[8:18]), tzinfo=timezone.utc)
        self.grids = []
        self.model = modelstr
        self.source = source
        self.trajectorycount = 0
        self.trajectories: list[Algorithm_HYSPLIT_Trajectory] = list()
        self.verticalmotion = None

    def addTrajectories(self, othert: list[Algorithm_HYSPLIT_Trajectory]) -> None:
        '''
        Add trajectories
        '''
        for t in othert:
            if not t in self.trajectories:
                self.trajectories.append(t)
                self.trajectorycount += 1

    def chart_time_slices(self, 
                          backgroundorbit: str,
                          backgroundprocessor: str, 
                          filter_start_start: datetime, 
                          filter_start_end: datetime, 
                          filter_end_start: datetime, 
                          filter_end_end: datetime, 
                          slice_date: 
                          datetime, fileimage:str) -> None:
        '''
        Chart emission/background air parcels

        Parameters
        ----------
        backgroundorbit: str
            TROPOMI orbit to chart under image
        filter_start_start: datetime
            Filter start for trajectory start
        filter_start_end: datetime
            Filter end for trajectory start
        filter_end_start: datetime
            Filter start for trajectory end
        filter_end_end: datetime
            Filter end for trajectory end
        slice_date: datetime
            Time of slice
        fileimage:str
            File to put output to. If None just display the image.
        '''
        geometry = Geometry.CONTAINS

        # Determine extent
        tbox = self.trajectories_extent(slice_date)
        [lonmin, latmin, lonmax, latmax] = tbox.bounds
        lonmin -= 0.1
        lonmax += 0.1
        latmin -= 0.1
        latmax += 0.1
        
        # Make it sligthly bigger
        tbox1 = box(lonmin, latmin, lonmax, latmax)

        # Determine proprtions of the image
        lonwidth = lonmax - lonmin
        latheight = latmax - latmin
        mll = lonwidth
        if mll < latheight: mll = latheight
        fracy = int(10 * latheight/mll) + 0.1
        if fracy < 4: fracy = 4

        geodetic = ccrs.Geodetic()
        plateCarree = ccrs.PlateCarree()
        fig = plt.figure(figsize=(10, fracy))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extents=(lonmin, lonmax, latmin, latmax), crs=plateCarree)
        ax.coastlines()

        if not backgroundorbit is None:
            tropomi: TROPOMI = TROPOMI_for_orbit(backgroundorbit, backgroundprocessor)
            [_, _, _, _, _, _, lons, lats, scanfiltered] = tropomi.narrow_to_domain(geometry, tbox1)

            v_ave = int(average(scanfiltered))
            v_max = v_ave + 25
            v_min = v_ave - 25
            cs = ax.pcolormesh(lons, lats, scanfiltered, shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.RdYlBu_r, transform=ccrs.PlateCarree())
            cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
            cbar.ax.set_xlabel("CH4 ppb")

        if tbox1.contains(self.source.xy):
            self.source.plot_source(ax)

        for t in self.trajectories:
            for (j,p) in t.points.items():
                if p.date == slice_date: 
                    p.draw_point(plt, geodetic, plateCarree)
                    if not logger is None: logger.info(p)

        title = "HYSPLIT Model: {}. Orbit: {} Processor: {}\nAir parcels at {}Z.".format(
            self.model,
            backgroundorbit, 
            tropomi.processor_version,
            slice_date.strftime("%Y%m%d %H"),
            )

        plt.title(title)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=2, color='gray', alpha=0.5, linestyle='--')

        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'} 
        

        if not v_ave is None: 
            fig.text(0.1, 0.08, "Box average: {} (ppb).".format(v_ave), horizontalalignment="left", wrap=False) 

        text5 = "Trajectories filtered start: [{}Z - {}Z] end: [{}Z - {}Z].".format(
            filter_start_start.strftime("%Y%m%d %H"), 
            filter_start_end.strftime("%Y%m%d %H"),
            filter_end_start.strftime("%Y%m%d %H"), 
            filter_end_end.strftime("%Y%m%d %H")
            )
        fig.text(0.1, 0.05, text5, horizontalalignment="left", wrap=False) 
        fig.subplots_adjust(bottom=0.2)

        if not fileimage is None:
            plt.savefig(fileimage)
            print("Chart generated: {}".format(fileimage))
            if not logger is None: logger.info("Chart generated: {}".format(fileimage))
        else:
            plt.show()

        plt.close("all")           

    def chart_time_series(self, geometry: Geometry, tropomiorbit: str, processor: str, fileimage: str):
        '''
        Chart containing HYSPLIT GDAS trajectories for a source on the top of TROPOMI image

        Parameters
        ----------
        geometry: Geometry
            Geometry used by TROPOMI
        tropomiorbi: str
            TROPOMI image to be overlayed with trajectories
        fileimage: str
            Path of the file to write the image to
        '''

        if 0 == len(self.trajectories):
            print("No HYSPLIT trajectories found")
            return

        # Determine extent
        latmax = -91
        latmin = 91
        lonmax = -181
        lonmin = 181

        for t in self.trajectories:
            for p in t.points.values():
                if p.longitude < lonmin: lonmin = p.longitude
                if lonmax < p.longitude : lonmax = p.longitude
                if p.latitude < latmin: latmin = p.latitude
                if latmax < p.latitude : latmax = p.latitude

        latmax += 0.1
        latmin -= 0.1
        lonmax += 0.1
        lonmin -= 0.1

        lonwidth = lonmax - lonmin
        latheight = latmax - latmin
        mll = lonwidth
        if mll < latheight: mll = latheight
        fracy = int(10 * latheight/mll)
        if fracy < 5: fracy = 5
        config = Config()
        tropomi: TROPOMI = TROPOMI_for_orbit(tropomiorbit, processor)
        [_, _, _, _, _, _, lons, lats, scanfiltered] = tropomi.narrow_to_domain(geometry, box(lonmin, latmin, lonmax, latmax))

        v_max = max(scanfiltered)
        v_min = min(scanfiltered)

        # print(v_min, v_max)
        # Overwrite
        v_min = 1800
        v_max = 1830
        geodetic = ccrs.Geodetic()
        plateCarree = ccrs.PlateCarree()
        fig = plt.figure(figsize=(10,fracy))

        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extents=(lonmin, lonmax, latmin, latmax), crs=plateCarree)
        ax.coastlines()

        cs = ax.pcolormesh(lons, lats, scanfiltered, shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.RdYlBu_r, transform=ccrs.PlateCarree())
        cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
        cbar.ax.set_xlabel("CH4 ppb")

        for t in self.trajectories:
            t.draw_trajectory(plt, geodetic, plateCarree, self.endtime)

        self.source.plot_source(ax)

        title = "HYSPLIT Model: {}. Orbit: {} Processor:{}\nBackward trajectory 2019091404 - 2019091423. Forward trajectory 2019091423 - 2019091504".format(self.model, tropomi.orbit, tropomi.processor_version)
        plt.title(title)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')

        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'} 

        if not fileimage is None:
            plt.savefig(fileimage)
            print("Chart generated: {}".format(fileimage))
        else:
            plt.show()

        plt.close("all")

    def chart_endpoints_count_per_pixel(self, geometry: Geometry, tropomiorbit: str,  processor:str, fileimage: str):
        '''
        This plots pixels that contains ends of trajectories and colours these by count
        '''

        endpoints = []
        for t in self.trajectories:
            for k,p in t.points.items():
                if (k == self.endtime): 
                    endpoints.append(Point(p.longitude, p.latitude))

        # overlay end of trajectories at 

        if 0 == len(self.trajectories):
            print("No HYSPLIT trajectory found")
            return

        # Determine extent
        latmax = -91
        latmin = 91
        lonmax = -181
        lonmin = 181
        for t in self.trajectories:
            for p in t.points.values():
                if p.longitude < lonmin: lonmin = p.longitude
                if lonmax < p.longitude : lonmax = p.longitude
                if p.latitude < latmin: latmin = p.latitude
                if latmax < p.latitude : latmax = p.latitude

        latmax += 0.1
        latmin -= 0.1
        lonmax += 0.1
        lonmin -= 0.1
        tropomi = TROPOMI_for_orbit(tropomiorbit, processor)
        [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, lons, lats, _]  = tropomi.narrow_to_domain(geometry, box(lonmin, latmin, lonmax, latmax))
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        endoftrajectoriescount = zeros((scans, pixels), dtype=int8)
        for s in range(scans):
            for p in range(pixels):
                ul = (lonbox[s,p,3], latbox[s,p,3])
                bl = (lonbox[s,p,0], latbox[s,p,0])
                br = (lonbox[s,p,1], latbox[s,p,1])
                ur = (lonbox[s,p,2], latbox[s,p,2])
                pixelbox = Polygon([ul,bl,br,ur])
                prepare(pixelbox)
                for a in endpoints:
                    if pixelbox.contains(a):
                        endoftrajectoriescount[s,p] += 1

        geodetic = ccrs.Geodetic()
        plateCarree = ccrs.PlateCarree()
        cmap = (ListedColormap(['cyan', 'blue', 'green']).with_extremes(under='white', over='red'))
        bounds=[1,2,3]
        norm = BoundaryNorm(bounds, cmap.N, extend="max")
        red_patch = mpatches.Patch(color="red", label="4 or more")
        green_patch = mpatches.Patch(color="green", label="3")
        blue_patch = mpatches.Patch(color="blue", label="2")
        cyan_patch = mpatches.Patch(color="cyan", label="1")

        # sm = ScalarMappable(norm = norm, cmap = cmap)
        fig = plt.figure(figsize=(10,10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extents=(lonmin, lonmax, latmin, latmax), crs=plateCarree)
        ax.coastlines()

        cs = ax.pcolormesh(lons, lats, endoftrajectoriescount, shading="nearest", norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
        ax.legend(handles=[red_patch, green_patch, blue_patch, cyan_patch])
        # cbar = fig.colorbar(ScalarMappable(norm=norm, cmap = cmap), ax=ax, orientation="horizontal", pad=0.075)

        self.source.plot_source(ax)

        for t in self.trajectories:
            t.draw_trajectory(plt, geodetic, plateCarree, self.endtime)

        title = "End of trajectoris count\n HYSPLIT Model {}\nOrbit {}".format(self.model, tropomiorbit)
        plt.title(title)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=2, color='gray', alpha=0.5, linestyle='--')

        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'} 

        plt.savefig(fileimage)
        plt.close("all")
        print("Generated chart {}".format(fileimage))
        pass

    def points_extent(self, polygon: Polygon, image_date: datetime) -> list[Algorithm_HYSPLIT_TrajectoryPoint]:
        '''
        Get points contained within polygon at specified time

        Parameters
        ----------
        polygon: Polygon 
            Emissions initial date
        image_date: datetime
            Time of trajectory points to analyse
        
        Returns
        -------
        list: list[Algorithm_HYSPLIT_TrajectoryPoint]
            List of points contained within the box at specified time
        '''
        results: list[Algorithm_HYSPLIT_TrajectoryPoint] = []
        prepare(polygon)
        if not logger is None: logger.info("Image date: {}".format(image_date))
        for t in self.trajectories:
            for p in t.points.values():
                if not polygon.contains(Point(p.longitude, p.latitude)): continue
                if not image_date is None and image_date != p.date: continue
                if not logger is None: logger.info("Point within domain: {}".format(p))
                results.append(p)

        if not logger is None: logger.info("{} points in the box".format(len(results)))
        return results 

    def trajectories_extent(self, slice_date: datetime) -> Polygon:
        '''
        Determine extent covering trajectories

        Parameters
        ----------
        slice_date: datetime
            Time of trajectory points to analyse
        
        Returns
        -------
        box: Polygon
            Polygon containing points at the time of image
        '''
        latmax = -91
        latmin = 91
        lonmax = -181
        lonmin = 181
        for t in self.trajectories:
            for p in t.points.values():
                if slice_date is None or slice_date != p.date: continue
                if p.longitude < lonmin: lonmin = p.longitude
                if lonmax < p.longitude : lonmax = p.longitude
                if p.latitude < latmin: latmin = p.latitude
                if latmax < p.latitude : latmax = p.latitude

        latmax += 0.1
        latmin -= 0.1
        lonmax += 0.1
        lonmin -= 0.1

        if not logger is None: logger.info("Trajectories box: [{}, {}] [{}, {}]".format(lonmin, latmin, lonmax, latmax))
        return box(lonmin, latmin, lonmax, latmax) 

def Algorithm_HYSPLIT_from_trajectory(source:Source, modelstr:str, datestr:str, filepath: str) -> Algorithm_HYSPLIT:
        '''
        Factory method constructing an Algorithm_HYSPLIT object from a file containing trajectories.
        Pattern for names of trajectory files. HYSPLIT_GFSQ_F_2019091405_2019091504.txt
        The file(s) should be located in the Source, Data,  Orbit, HYSPLIT folder
        1. HYSPLIT_
        2. model - CDC, GDAS or GFSQ
        3. direction - backwards B of forwards F
        5. Trajectory start time - YYYYMMDDHH
        5. Trajectory end time - YYYYMMDDHH
        6. Extension - .txt

        Parameters
        ----------
        source: Source
            Source of emissions
        modelstr: str
            Type of data used by HYSPLIT "CDC", "GDAS", "GFSQ"
        datestr: str
        filepath: str
            Path to the trajectory file

        '''
        h = Algorithm_HYSPLIT(source, modelstr, datestr)
        numberofgrids = 0
        record_1_start = 1
        record_2_start = 2
        record_3_start = 0
        record_4_start = 0
        record_5_start = 0
        record_6_start = 0
        linenumber = 1
        trajectorycounter = 1
        if not path.exists(filepath):
            print("Error: Unable to open file {}. Please run HYSPLIT and create trajectory file first.".format(filepath))
            return None
        
        with open(filepath, 'r', encoding='UTF-8') as f:
            lines = f.readlines()
            f.close()

        if 0 == len(lines):
            print("Error: File {} is empty".format(filepath))
            return None

        for line in lines:
            s = line.strip()
            arr = s.split()

            # Record 1
            if record_1_start == linenumber:
                numberofgrids = int(arr[0])
                record_3_start = record_2_start + numberofgrids
                record_4_start = record_3_start + 1
                linenumber += 1

            # Loop Records #2 
            elif record_2_start <= linenumber and linenumber < record_3_start:
                h.grids.append(arr)
                linenumber += 1

            # Record #3
            elif record_3_start == linenumber:
                h.trajectorycount = int(arr[0])
                h.direction = arr[1]
                h.verticalmotion = arr[2]
                record_5_start = record_4_start + h.trajectorycount
                record_6_start = record_5_start + 1
                linenumber += 1

            # Loop Record #4 ==> number of different trajectories in file
            elif record_4_start <= linenumber and linenumber < record_5_start:
                t_start = datetime(int(arr[0])+2000, int(arr[1]), int(arr[2]), int(arr[3]), tzinfo=timezone.utc)
                t_latitude = float(arr[4])
                t_longitude = float(arr[5])
                t_height = float(arr[6])
                h.trajectories.append(Algorithm_HYSPLIT_Trajectory(t_latitude, t_longitude, t_height))
                trajectorycounter += 1
                linenumber += 1

            # Record #5
            elif record_5_start == linenumber:
                numberofdiagnostics = int(arr[0])
                h.diagnostics = arr[1:]
                linenumber += 1

            # Loop Record #6 ==> through end of all endpoints
            elif record_6_start <= linenumber:
                t = int(arr[0]) - 1
                gridid = int(arr[1])
                date = datetime(2000 + int(arr[2]), int(arr[3]), int(arr[4]), int(arr[5]), int(arr[6]), tzinfo=timezone.utc)
                forecasthour = int(arr[7])
                ageoftrajectory = float(arr[8])
                latitude = float(arr[9])
                longitude = float(arr[10])
                height = float(arr[11])
                diagnostics: list[float] = []
                for s in arr[12:]: diagnostics.append(float(s))
                hpoint = Algorithm_HYSPLIT_TrajectoryPoint(gridid, date, forecasthour, ageoftrajectory, latitude, longitude, height, diagnostics)
                h.trajectories[t].points[date] = hpoint
                linenumber += 1
        
        return h

def _handler_chart_time_series(args: Namespace):
    '''
    Driver for charting function of HYSPLIT algorithm
    '''
    Config.create_log("Algorithm_HYSPLIT_{}_{}.log".format(args.source, args.orbit))

    if args.source == "HailCreek":
        source = create_source("HailCreek")

    geometry = Geometry.CENTER

    ah = Algorithm_HYSPLIT(source, args.model, args.enddate)
    t_backward = path.join(Config.HYSPLIT_folder, args.source, "HYSPLIT_{}_B_{}_{}.txt".format(args.model, args.startdate, args.date_at_mine))
    t_forward = path.join(Config.HYSPLIT_folder, args.source, "HYSPLIT_{}_F_{}_{}.txt".format(args.model, args.date_at_mine, args.enddate))
    if not logger is None: logger.info("Processing file {}".format(t_backward))
    if path.exist(t_backward): ah.addTrajectories(Algorithm_HYSPLIT_from_trajectory(source, args.model, args.date_at_mine, t_backward).trajectories)
    if path.exist(t_forward): ah.addTrajectories(Algorithm_HYSPLIT_from_trajectory(source, args.model, args.date_at_mine, t_forward).trajectories)

    ah.chart_time_series(geometry, args.orbit, args.processor, None)

def _handler_chart_time_slice(args: Namespace):
    '''
    Chart containing HYSPLIT trajectories for a source on the top of TROPOMI image
    '''
    Config.create_log("Algorithm_HYSPLIT_{}_{}.log".format(args.source, args.orbit))

    sourcepath = path.join(Config.HYSPLIT_folder, args.source)
    filenames = listdir(sourcepath)
    pattern = "HYSPLIT_{}_{}_".format(args.model, args.direction)
    source = create_source(args.source)
    filter_start_start = datetime(int(args.filter_start_start[0:4]), int(args.filter_start_start[4:6]), int(args.filter_start_start[6:8]), int(args.filter_start_start[8:10]), tzinfo=timezone.utc)
    filter_start_end = datetime(int(args.filter_start_end[0:4]), int(args.filter_start_end[4:6]), int(args.filter_start_end[6:8]), int(args.filter_start_end[8:10]), tzinfo=timezone.utc)
    filter_end_start = datetime(int(args.filter_end_start[0:4]), int(args.filter_end_start[4:6]), int(args.filter_end_start[6:8]), int(args.filter_end_start[8:10]), tzinfo=timezone.utc)
    filter_end_end = datetime(int(args.filter_end_end[0:4]), int(args.filter_end_end[4:6]), int(args.filter_end_end[6:8]), int(args.filter_end_end[8:10]), tzinfo=timezone.utc)
    slice_datetime = datetime(int(args.date[0:4]), int(args.date[4:6]), int(args.date[6:8]), int(args.date[8:10]), tzinfo=timezone.utc)

    ah = Algorithm_HYSPLIT(source, args.model, slice_datetime)
    for f in filenames:
        if not f.startswith(pattern): continue
        if not logger is None: logger.info("Processing file {}".format(f))
        trajectories = Algorithm_HYSPLIT_from_trajectory(source, "GFSQ", args.date, path.join(sourcepath, f)).trajectories
        filtered = []
        for t in trajectories:
            e_dt = t.end_datetime()
            s_dt = t.start_datetime()
            if not logger is None: 
                logger.info("Start:{} End:{} Chart: {}".format(s_dt.strftime("%Y%m%d%H"), e_dt.strftime("%Y%m%d%H"), slice_datetime.strftime("%Y%m%d%H")))
            if filter_start_start <= s_dt and s_dt <= filter_start_end and filter_end_start <= e_dt and e_dt <= filter_end_end:
                filtered.append(t)
                if not logger is None: logger.info("Adding trajectory from file: {}".format(f))
        ah.addTrajectories(filtered)

    if 0 == len(ah.trajectories):
        if not logger is None: logger.error("No HYSPLIT trajectory files found")
        print("No HYSPLIT trajectory files found")
        return

    ah.chart_time_slices(args.backgroundorbit, args.backgroundprocessor,filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, None)

def _handler_clean(args):
    '''
    Handler for clean menu
    '''
    logname = "Algorithm_HYSPLIT_{}_{}.log".format(args.source, args.orbit)
    filename_log = path.join("Log", logname)

    if path.exists(filename_log): remove(filename_log)
    return

def _handler_create_files(args):
    '''
    Handler for creating files for trajectories for time between the orbit and previous orbit
    Previous orbit must be less then 26 hours
    '''
    assets_folder = path.join(Config.HYSPLIT_folder, args.source)
    source = create_source(args.source)
    tropomi_emissions = TROPOMI_for_orbit(args.orbit, args.processor)
    tropomi_background = TROPOMI_for_orbit(args.backgroundorbit, args.backgroundprocessor)
    [_, _, _, emission_date] = tropomi_emissions.get_pixel_for_source(source.xy)
    [_, _, _, background_date] = tropomi_background.get_pixel_for_source(source.xy)

    if background_date + timedelta(hours=26) < emission_date or emission_date < background_date + timedelta(hours=1):
        m = "Error: Previous orbit must be more then an hour before and less than 26 hours before {} {}".format(background_date.strftime("%Y%m%d%H"), emission_date.strftime("%Y%m%d%H"))
        print(m)
        if not logger is None: logger.error(m)
        return
    
    if not path.exists(assets_folder): Config.create_directory_structure_HYSPLIT(args.source)
    pattern = "HYSPLIT_{}_{}_{}_{}.txt"

    dt = emission_date - timedelta(hours = 1)
    while background_date < dt:
        filename = path.join(assets_folder, pattern.format(args.model, "B", background_date.strftime("%Y%m%d%H"), dt.strftime("%Y%m%d%H")))
        if not path.exists(filename):
            with open(filename, "at") as f:
                f.flush()
                f.close
                print(filename)
        dt -= timedelta(hours = 1)

    dt = background_date
    while dt < emission_date:
        filename = path.join(assets_folder, pattern.format(args.model, "F", dt.strftime("%Y%m%d%H"), emission_date.strftime("%Y%m%d%H")))
        if not path.exists(filename):
            with open(filename, "at") as f:
                f.flush()
                f.close
                print(filename)
        dt += timedelta(hours = 1)
    pass

if __name__ == '__main__':
    parser = ArgumentParser(prog="Algorithm_HYSPLIT", description="Manage HYSPLIT trajectories")
    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_chart= subparsers.add_parser("chart", help="Chart hysplit")
    subparser_charts = parser_chart.add_subparsers(title="subcommands")

    parser_chart_trajectory = subparser_charts.add_parser("series", help="Chart air parcel time series")
    parser_chart_trajectory.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions", )
    parser_chart_trajectory.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 09956 used for location of HYSPLIT file")
    parser_chart_trajectory.add_argument("processor", type=str, help="Processor for orbit 020400.")
    parser_chart_trajectory.add_argument("model", type=str, choices=["CDC", "GDAS", "GFSQ"])
    parser_chart_trajectory.add_argument("startdate", type=str, help="Start date inclusive in the format YYYYMMDDHH e.g. 2019091504")
    parser_chart_trajectory.add_argument("date_at_mine", type=str, help="Date at mine inclusive in the format YYYYMMDDHH e.g. 2019091504")
    parser_chart_trajectory.add_argument("enddate", type=str, help="Last date inclusive in the format YYYYMMDDHH e.g. 2019091504")
    parser_chart_trajectory.set_defaults(func=_handler_chart_time_series)

    parser_chart_puffs = subparser_charts.add_parser("slice", help="Chart air parcels time slice")
    parser_chart_puffs.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions", )
    parser_chart_puffs.add_argument("orbit", type=str, help="Orbit of emissions as a five digit string e.g. 09956")
    parser_chart_puffs.add_argument("processor", type=str, help="Processor for orbit 020400.")
    parser_chart_puffs.add_argument("backgroundorbit", type=str, help="Orbit of background as a five digit string e.g. 09942")
    parser_chart_puffs.add_argument("backgroundprocessor", type=str, help="Processor for background orbit 020400.")
    parser_chart_puffs.add_argument("model", type=str, choices=["CDC", "GDAS", "GFSQ"])
    parser_chart_puffs.add_argument("direction", type=str, choices=["F", "B"], help="Forward or backward trajectories")
    parser_chart_puffs.add_argument("filter_start_start", type=str, help="Trajectory start filter start date inclusive in the format YYYYMMDDHH e.g. 2019091504")
    parser_chart_puffs.add_argument("filter_start_end", type=str, help="Trajectory start filter end date inclusive in the format YYYYMMDDHH e.g. 2019091504")
    parser_chart_puffs.add_argument("filter_end_start", type=str, help="Trajectory end filter start date inclusive in the format YYYYMMDDHH e.g. 2019091504")
    parser_chart_puffs.add_argument("filter_end_end", type=str, help="Trajectory end filter end date inclusive in the format YYYYMMDDHH e.g. 2019091504")
    parser_chart_puffs.add_argument("date", type=str, help="Date of time slice in the format YYYYMMDDHH e.g. 2019091504")
    parser_chart_puffs.set_defaults(func=_handler_chart_time_slice)

    parser_clean = subparsers.add_parser("clean", help="Clean outputs for specified source and orbit.")
    parser_clean.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_clean.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 09956 used for location of HYSPLIT file")
    parser_clean.add_argument("model", type=str, choices=["CDC", "GDAS", "GFSQ"])
    parser_clean.set_defaults(func=_handler_clean)

    parser_create_files = subparsers.add_parser("create", help="Create files for hysplit trajectories")
    parser_create_files.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions", )
    parser_create_files.add_argument("orbit", type=str, help="Orbit of emissions as a five digit string e.g. 09956. ")
    parser_create_files.add_argument("processor", type=str, help="Processor for orbit 020400.")
    parser_create_files.add_argument("backgroundorbit", type=str, help="Orbit of background as a five digit string e.g. 09942. Used as a background image.")
    parser_create_files.add_argument("backgroundprocessor", type=str, help="Processor for background orbit 020400.")
    parser_create_files.add_argument("model", type=str, choices=["CDC", "GDAS", "GFSQ"])
    parser_create_files.set_defaults(func=_handler_create_files)

    args = parser.parse_args()
    
    if ("orbit" in args):
        if not(len(args.orbit) == 5) or not(args.orbit.isdigit()):
            print ("Orbit must be a five digit number: {}".format(args.orbit))
            exit()

    if ("backgroundorbit" in args):
        if not(len(args.backgroundorbit) == 5) or not(args.backgroundorbit.isdigit()):
            print ("Background orbit must be a five digit number: {}".format(args.backgroundorbit))
            exit()

    if ("processor" in args):
        if not(len(args.processor) == 6) or not(args.processor.isdigit()):
            print ("Processor must be a 6 digit number: {}".format(args.processor))
            exit()

    args.func(args)
