from argparse import ArgumentParser, Namespace
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy import geodesic
from Chemistry import convert_column_ppb_with_water
from Config import Config
from datetime import datetime
from enum import Enum
from Geometry import Geometry
from Mask import Mask
import logging
from math import cos, pi, sqrt
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from Meteorology import Meteorology
from numpy import average, ma, nan
from os import listdir, path, remove
from pickle import dump, load
from shapely import LineString,prepare, Point, Polygon, box
from Source import create_source, Source
from trajectory import drawTrajectory
from TROPOMI import TROPOMI, TROPOMI_for_orbit

logger = logging.getLogger(__name__)

class EnhancementLength:
    '''
    Total enhancement for length adjustment
    '''
    def __init__(self, length: float, enhancement: float, increment: float):
        self.length = length
        self.enhancement = enhancement
        self.increment = increment

    def __str__(self):
        if self.increment is None: return "Length: {}, Enhancement: {}, Increment: Not defined".format(self.length, self.enhancement)
        return "Length: {}, Enhancement: {}, Increment: {}".format(self.length, self.enhancement, self.increment)

class EnhancementWidth:
    '''
    Total enhancement for width adjustment
    '''
    def __init__(self, half_width: float, enhancement: float, increment: float):
        self.half_width = half_width
        self.enhancement = enhancement
        self.increment = increment

    def __str__(self):
        if self.increment is None: return "Half-width: {}, Enhancement: {}, Increment: Not defined".format(self.half_width, self.enhancement)
        return "Half-width: {}, Enhancement: {}, Increment: {}".format(self.half_width, self.enhancement, self.increment)

class Rotation:
    '''
    Average enhancement for rotation adjustment
    '''
    def __init__(self, rotation: float, azimuth: float, averageenhancement: float):
        self.rotation = rotation
        self.azimuth = azimuth
        self.averageenhancement = averageenhancement

    def __str__(self):
        return "Rotation: {}, Azimuth: {}, Average Enhancement: {}".format(self.rotation, self.azimuth, self.averageenhancement)
    
class Status(Enum):
    '''
    Algoritm status as a bit map
    '''
    PROCESSING = 0
    SUCCESS = 1
    FAILURE_TRANSECTS_BELOW_3 = 2 # There must be at least 3 valid transects
    FAILURE_WIND_BELOW2 = 4 # Wind must be greater then 2 m/s
    FAILURE_NODATA_DOWNWINDBOX = 8 # Downwind box contains no data
    FAILURE_NODATA_INDOMAIN = 16 # Domain contains no data
    FAILURE_BACKGROUND_VALUES_HIGH = 32 # Background contains no data
    FAILURE_TRANSECTS_PIXEL_COUNT = 64 # Number of valid pixels intersecting transects below threshold
    FAILURE_ROTATION_0_NO_PIXELS = 128 # Rotation 0 does not contain any pixels

class Transect:
    '''
    From Sadavarte 2021: We define 15 equally spaced transects between the source and the end of the rectangular mask for calculating the source rates.
    We ignore the first three transects due to their close proximity to the source, where XCH4 may be underestimated due to partial
    pixel enhancement. To avoid underestimation of emissions due to incomplete sampling of the plume by a transect due to 
    missing pixels, we only consider transects that have more than 75% overlap with TROPOMI pixels. With this requirement, we 
    only calculate the source rate from plumes with at least three or more transects
    '''
    def __init__(self, start_lon: float, start_lat: float, end_lon: float, end_lat: float):
        self.enhkg: float = 0.0
        self.len_mask: float = 0.0
        self.len_values: float = 0.0
        self.len_values_positive: float = 0
        self.pix_mask: int = 0
        self.pix_values: int = 0
        self.pix_values_positive: int = 0
        self.line: LineString = LineString([[start_lon, start_lat],[end_lon, end_lat]])

    def __repr__(self):
        return "Transect({},{},{},{})".format(self.line.coords[0][0], self.line.coords[0][1], self.line.coords[1][0], self.line.coords[1][1])
        
    def __str__(self) -> str:
        str0 = "Valid: {} ".format(self.isvalid())
        str1 = "Enh: {0: .3f}(kg/m) Length: [masked: {1:.3f}(m) valid: {2:.3f}(m) positive: {3:.3f}(m)] ". format(self.enhkg, self.len_mask, self.len_values, self.len_values_positive)
        str2 = "Count: [masked: {} valid: {} positive: {}] ".format(self.pix_mask, self.pix_values, self.pix_values_positive)
        str3 = "Line: [{:.3f},{:.3f}] - [{:.3f},{:.3f}]".format(self.line.coords.xy[0][0],self.line.coords.xy[1][0],self.line.coords.xy[0][1],self.line.coords.xy[1][1])
        return str0 + str1 + str2 + str3
        
    def isvalid(self) -> bool:
        '''
        From Sadavarte2021: To avoid underestimation of emissions due to incomplete sampling of the plume by a transect due to missing pixels, 
        we only consider transects that have more than 75% overlap with TROPOMI pixels.
        '''
        if (self.pix_values < 1): return False
        # Added sanity condition that there must be at least one positive pixel. This is not in the orignal paper
        if (self.pix_values_positive < 1): return False
        # Interpret validity in terms of pixel count
        if (self.pix_values < self.pix_mask * 3 ): return False
        # Interpret validity in terms of length
        if (self.len_values < self.len_mask * 3 ): return False
        return True

class Algorithm_CSF:
    '''
    Implementation of Cross Sectional Flux (CSF) method CH4 emission estimation algorithm as described in: 
    Methane Emissions from Superemitting Coal Mines in Australia Quantified Using TROPOMI Satellite Observations
    Pankaj Sadavarte, Sudhanshu Pandey, Joannes D. Maasakkers, Alba Lorente, Tobias Borsdorff, Hugo Denier van der Gon, Sander Houweling, and Ilse Aben
    Environmental Science & Technology 2021 55 (24), 16573-16580
    DOI: 10.1021/acs.est.1c03976

    These calculations require the following inputs:
    1. Location of mine
    2. Pressure boundarly layer height averaged wind at the time of TROPOMI image at location of the mine
    3. Surface pressure a the time of TROPOMI image covering the plume
    '''
    
    degree = 110000.
    degree_helf = degree / 2.
    degree_quarter = degree / 4.
    degree_tenth = degree / 10.
    degree_2tenth = 2. * degree_tenth
    degree_4tenth = 4. * degree_tenth
    degree_twentieth = degree / 20.
    downwind_box_length_initial = degree_4tenth
    downwind_box_width_initial = degree_2tenth
    
    g = geodesic.Geodesic()

    def __init__(self, config: Config, source: Source, wind: list[float], tropomi: TROPOMI):
        '''
        Initialize algorithm

        Parameters
        ----------
        config: Config
            Configuration for algorithm run
        source: Source
            Source of emissions
        wind: float[u,v]
            Eastward component of pressure averaged boundary layer wind
            Northward component of pressure averaged boundary layer wind
        tropomi: TROPOMI
            TROPOMI file containing CH4 scan, pixel center and pixel box coordinates
        '''
        if (not(logger is None)): logger.info("Source: lon:{} lat: {}".format(source.longitude, source.latitude))
        if (not(logger is None)): logger.info(config)

        self.config = config

        self.background: float = None
        self.background_enhancement: float = None

        self.domain_median: float = None

        self.downwind_box: Polygon = None
        self.downwind_box_azimuth_adjusted: float = None
        # As the algoritm calls for extension of the downwind box in the direction of wind we used ocenaographic convention here
        self.downwind_box_azimuth_initial: float = Meteorology.calculate_azimuth_degree_oceanography(wind[0], wind[1])
        self.downwind_box_count = 0
        self.downwind_box_endpoint: Point = source.xy
        self.downwind_box_half_width: float = Algorithm_CSF.degree_tenth
        self.downwind_box_initial: Polygon = None
        self.downwind_box_length: float = Algorithm_CSF.degree_4tenth
        self.downwind_box_minimumvalue: float = None
        self.enhancement_length: list[EnhancementLength] = []
        self.enhancement_width: list[EnhancementWidth] = []
        self.rotations: list[Rotation] = []
        self.shape_enhancement: float = None

        self.source = source
        
        self.transects: list[Transect] = None
        self.transects_valid_count = 0
        self.transects_valid_pixels_count = 0 # Number of valid pixels intersected by transects
        self.transects_valid_positive_pixels_count = 0 # Number of valid pixels with positive values intersected by transects

        self.tropomi: TROPOMI = tropomi
        self.tropomi_source_scan: int = None
        self.tropomi_source_pixel: int = None
        self.tropomi_source_date: datetime = None

        self.upwind_point_01: Point = None
        self.upwind_box = None
        self.upwind_box_background_average = None
        self.upwind_box_background_count = None

        self.wind = wind
        self.u: float = wind[0]
        self.v: float = wind[1]

        # Calculations
        self.q_valid_hour: float = None # Emission rate t/hour
        self.q_valid_year: float = None # Emission rate per year Gg/year
        self.status: int = Status.PROCESSING.value # Status of processing

        # Initialize values
        [self.tropomi_source_scan, self.tropomi_source_pixel, self.tropomi_source_date, _] = tropomi.get_pixel_for_source(self.source.xy)

        if (self.tropomi_source_scan is None) or (self.tropomi_source_pixel is None) or (self.tropomi_source_date is None): 
            message = "Unable to locate source within tropomi scan for orbit {}. Program terminates.".format(self.tropomi.orbit)
            if not (logger is None): logger.error(message)
            print(message)
            exit()      

        # Create folder structure just in case it does not exists
        Config.create_directory_structure_CSF(source.case_name)

        self.picklename = Algorithm_CSF.get_picklename(self.config, source.case_name, tropomi.orbit, tropomi.processor_version)

        # Calculate initial downwind box
        p1 = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=self.downwind_box_azimuth_initial, distances=self.downwind_box_length)
        pul = Algorithm_CSF.g.direct(points=[p1[0,0], p1[0,1]], azimuths=self.downwind_box_azimuth_initial - 90, distances=Algorithm_CSF.degree_tenth)
        pbl = Algorithm_CSF.g.direct(points=[p1[0,0], p1[0,1]], azimuths=self.downwind_box_azimuth_initial + 90, distances=Algorithm_CSF.degree_tenth)
        pur = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=self.downwind_box_azimuth_initial-90, distances=Algorithm_CSF.degree_tenth)
        pbr = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=self.downwind_box_azimuth_initial+90, distances=Algorithm_CSF.degree_tenth)
        self.downwind_box_initial = Polygon([(pul[0,0], pul[0,1]), (pbl[0,0], pbl[0,1]), (pbr[0,0], pbr[0,1]), (pur[0,0], pur[0,1])])

        if (not(logger is None)): logger.info("Downwind box initial: [[{}, {}], [{}, {}], [{}, {}], [{}, {}]]".format(
            pul[0,0], pul[0,1], pbl[0,0], pbl[0,1], pbr[0,0], pbr[0,1], pur[0,0], pur[0,1]))

        # Calculate upwind box
        # As the algoritm calls for extension of the upwind box in the direction of upwind wind we used ocenaographic convention here
        azimuth_upwind = Meteorology.calculate_azimuth_degree_oceanography(-self.u, -self.v) 
        p1_u = Algorithm_CSF.g.direct(points=[self.source.xy.x, self.source.xy.y], azimuths=azimuth_upwind, distances=Algorithm_CSF.degree_tenth)
        pul_u = Algorithm_CSF.g.direct(points=[p1_u[0,0], p1_u[0,1]], azimuths=azimuth_upwind - 90, distances=Algorithm_CSF.degree_quarter)
        pbl_u = Algorithm_CSF.g.direct(points=[p1_u[0,0], p1_u[0,1]], azimuths=azimuth_upwind + 90, distances=Algorithm_CSF.degree_quarter)
        pur_u = Algorithm_CSF.g.direct(points=[pul_u[0,0], pul_u[0,1]], azimuths=azimuth_upwind, distances=Algorithm_CSF.degree_helf)
        pbr_u = Algorithm_CSF.g.direct(points=[pbl_u[0,0], pbl_u[0,1]], azimuths=azimuth_upwind, distances=Algorithm_CSF.degree_helf)
        self.upwind_point_01 = Point(p1_u[0,0], p1_u[0,1])
        self.upwind_box = Polygon([(pul_u[0,0], pul_u[0,1]), (pbl_u[0,0], pbl_u[0,1]), (pbr_u[0,0], pbr_u[0,1]), (pur_u[0,0], pur_u[0,1])])

        if (not(logger is None)): logger.info("Upwind box initial: [[{}, {}], [{}, {}], [{}, {}], [{}, {}]]".format(
            pul[0,0], pul[0,1], pbl[0,0], pbl[0,1], pbr[0,0], pbr[0,1], pur[0,0], pur[0,1]))

        # Save state into a file
        with open(self.picklename, mode="wb") as f:
            dump(self, f)
            f.close()

    def _calculate_background(self):
        '''
        From Sadavarte 2021: The methane enhancements for each pixel along the transects is defined relative to the 
        background XCH4 which is calculated as the average of 0.5° x 0.5° area centered at a distance of 
        0.1° upwind from the source. 

        From Sadavarte 2021: If the number of background observations is less than 20, we use the 
        median XCH4 of all pixels in the domain (20°-24°S, 146°-150°E) as background XCH4. To 
        account for other emissions in the downwind plume, we subtract the contributions from 
        surrounding coal mines

        On exit this function sets self.background is set 
        '''
        if (not(logger is None)): logger.info("")

        [minscan, minpixel, maxscan, maxpixel, lon_box, lat_box, lon_ctr, lat_ctr, scan_domain] = self.tropomi.narrow_to_domain(self.config.Algorithm_CSF_background_geometry, self.config.Algorithm_CSF_domain)
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        if (scan_domain is None):
            self.status = Status.FAILURE_NODATA_INDOMAIN.value
            with open(self.picklename, mode="wb") as f:
                dump(self, f)
                f.close()
            return()
        
        for s in range(scans):
            for p in range(pixels):
                ul = (lon_box[s,p,3], lat_box[s,p,3])
                bl = (lon_box[s,p,0], lat_box[s,p,0])
                br = (lon_box[s,p,1], lat_box[s,p,1])
                ur = (lon_box[s,p,2], lat_box[s,p,2])
                pixelbox = Polygon([ul,bl,br,ur])
                center = Point(lon_ctr[s,p], lat_ctr[s,p])
                if self.config.Algorithm_CSF_background_geometry == Geometry.CENTER and not(self.config.Algorithm_CSF_domain.contains(center)): scan_domain[s,p] = ma.masked
                if self.config.Algorithm_CSF_background_geometry == Geometry.CONTAINS and not(self.config.Algorithm_CSF_domain.contains(pixelbox)): scan_domain[s,p] = ma.masked
                if self.config.Algorithm_CSF_background_geometry == Geometry.INTERSECTS and not(self.config.Algorithm_CSF_domain.intersects(pixelbox)): scan_domain[s,p] = ma.masked

        self.domain_median = ma.median(scan_domain)
        if (not(logger is None)): logger.info("Domain median: {} ppb".format(self.domain_median))

        prepare(self.upwind_box)
        [minscan, minpixel, maxscan, maxpixel, lon_box, lat_box, lon_ctr, lat_ctr, scand_upwindbox] = self.tropomi.narrow_to_domain(self.config.Algorithm_CSF_background_geometry, self.upwind_box)
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        for s in range(scans):
            for p in range(pixels):
                ul = (lon_box[s,p,3], lat_box[s,p,3])
                bl = (lon_box[s,p,0], lat_box[s,p,0])
                br = (lon_box[s,p,1], lat_box[s,p,1])
                ur = (lon_box[s,p,2], lat_box[s,p,2])
                pixelbox = Polygon([ul,bl,br,ur])
                center = Point(lon_ctr[s,p], lat_ctr[s,p])
                if self.config.Algorithm_CSF_background_geometry == Geometry.CENTER and not(self.upwind_box.contains(center)): scand_upwindbox[s,p] = ma.masked
                if self.config.Algorithm_CSF_background_geometry == Geometry.CONTAINS and not(self.upwind_box.contains(pixelbox)): scand_upwindbox[s,p] = ma.masked
                if self.config.Algorithm_CSF_background_geometry == Geometry.INTERSECTS and not(self.upwind_box.intersects(pixelbox)): scand_upwindbox[s,p] = ma.masked

        self.upwind_box_background_average = ma.average(scand_upwindbox)
        self.upwind_box_background_count = ma.count(scand_upwindbox)

        if (self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count):
            self.background_enhancement = self.upwind_box_background_average
            if (not(logger is None)): logger.info("ALGORITHM: Background Enhancement: Upwind box average. Masked count: {} average: {}".format(self.upwind_box_background_count, self.upwind_box_background_average))
        else:
            if (not(logger is None)): logger.info("ALGORITHM: Background Enhancement: Domain median: {}".format(self.domain_median))
            self.background_enhancement = self.domain_median

        if (not(logger is None)): logger.info("Upwind box average : {} ppb".format(self.upwind_box_background_average))

        if self.config.Algorithm_CSF_transect_background == "algorithm":
            if (self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count):
                self.background = self.upwind_box_background_average
                if (not(logger is None)): logger.info("ALGORITHM: Upwind box average. Masked count: {} average: {}".format(self.upwind_box_background_count, self.upwind_box_background_average))
            else:
                if (not(logger is None)): logger.info("ALGORITHM: Domain median: {}".format(self.domain_median))
                self.background = self.domain_median

        if self.config.Algorithm_CSF_transect_background == "domain": 
                if (not(logger is None)): logger.info("OPTION: Domain median: {}".format(self.domain_median))
                self.background = self.domain_median

        if self.config.Algorithm_CSF_transect_background == "manual": 
                if (not(logger is None)): logger.info("OPTION: Manual average: {}".format(self.config.Algorithm_CSF_transect_background_value))
                self.background  = self.config.Algorithm_CSF_transect_background_value

        if self.config.Algorithm_CSF_transect_background == "upwind": 
                if (not(logger is None)): logger.info("OPTION: Upwind box average. Masked count: {} average: {}".format(self.upwind_box_background_count, self.upwind_box_background_average))
                self.background  = self.config.Algorithm_CSF_transect_background_value

        with open(self.picklename, mode="wb") as f:
            dump(self, f)
            f.close()

        return

    def _calculate_transects(self):
        '''
        Private interface to calculate_transects
        '''
        self.calculate_transects()

        with open(self.picklename, mode="wb") as f:
            dump(self, f)
            f.close()
        
        return

    def _calculate_downwindbox_directionadjustment(self):
        '''
        From Sadavarte 2021: We start with a smaller rectangular mask of dimension (length x breadth) 0.4 x 0.2° placed at the source in the
        downwind direction inferred from boundary layer average ERA5 meteorology to define the area containing the plume (Figure S1). 

        Next, we rotate this mask from -40 to +40° at 5° intervals around the inferred ERA5 wind direction such that the average
        XCH4 enhancement in the rectangular mask is maximal.
        
        On exit set self.pointDownindEnd and self.DownwindBox
        '''
        if not (self.status == Status.PROCESSING.value): return
        if (None == self.background): self._calculate_background()

        # As default value use azimuth = 0 
        azimuth = self.downwind_box_azimuth_initial
        p1 = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=azimuth, distances=self.downwind_box_length)
        pnt = Point(p1[0,0], p1[0,1])
        pul = Algorithm_CSF.g.direct(points=pnt.coords, azimuths=azimuth - 90, distances=Algorithm_CSF.degree_tenth)
        pbl = Algorithm_CSF.g.direct(points=pnt.coords, azimuths=azimuth + 90, distances=Algorithm_CSF.degree_tenth)
        pur = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=azimuth-90, distances=Algorithm_CSF.degree_tenth)
        pbr = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=azimuth+90, distances=Algorithm_CSF.degree_tenth)
        poly = Polygon([(pul[0,0], pul[0,1]), (pbl[0,0], pbl[0,1]), (pbr[0,0], pbr[0,1]), (pur[0,0], pur[0,1])])
        
        [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, lons, lats, scan] = self.tropomi.narrow_to_domain(self.config.Algorithm_CSF_downwindbox_geometry, poly)

        if (maxscan is None) or (minscan is None) or (maxpixel is None) or (minpixel is None): 
            if (not(logger is None)): logger.info("Rejected: rotation: 0 (deg). Scan does not contain pixels.")
            self.status = Status.FAILURE_ROTATION_0_NO_PIXELS.value

            with open(self.picklename, mode="wb") as f:
                dump(self, f)
                f.close()            
            return()

        scans = maxscan - minscan
        pixels = maxpixel - minpixel

        [self.downwind_box_count, self.shape_enhancement, old_average] = self._calculate_downwindbox_enhancement_statistics(poly, scans, pixels, lonbox, latbox, lons, lats, scan)

        self.downwind_box_endpoint = pnt
        self.downwind_box = poly
        self.downwind_box_azimuth_adjusted = azimuth
        if (not(logger is None)): logger.info("Accepted: default rotation: 0 (deg) average: {} azimuth {}".format(old_average , azimuth))
        self.rotations.append(Rotation(0, azimuth, old_average))
        for r in range(-40,45,5):
            if r == 0: continue

            azimuth = self.downwind_box_azimuth_initial + r
            p1 = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=azimuth, distances=self.downwind_box_length)
            pnt = Point(p1[0,0], p1[0,1])
            pul = Algorithm_CSF.g.direct(points=pnt.coords, azimuths=azimuth - 90, distances=Algorithm_CSF.degree_tenth)
            pbl = Algorithm_CSF.g.direct(points=pnt.coords, azimuths=azimuth + 90, distances=Algorithm_CSF.degree_tenth)
            pur = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=azimuth-90, distances=Algorithm_CSF.degree_tenth)
            pbr = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=azimuth+90, distances=Algorithm_CSF.degree_tenth)
            poly = Polygon([(pul[0,0], pul[0,1]), (pbl[0,0], pbl[0,1]), (pbr[0,0], pbr[0,1]), (pur[0,0], pur[0,1])])
            [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, lons, lats, scan]  = self.tropomi.narrow_to_domain(self.config.Algorithm_CSF_downwindbox_geometry, poly)

            if maxscan is None or minscan is None or maxpixel is None or minpixel is None: 
                if (not(logger is None)): logger.info("Rejected: rotation: {} Scan does not contain pixels".format(r))
                continue

            scans = maxscan - minscan
            pixels = maxpixel - minpixel


            [c, s, average] = self._calculate_downwindbox_enhancement_statistics(poly, scans, pixels, lonbox, latbox, lons, lats, scan)

            if ((old_average is None) or (old_average < average )):
                old_average = average
                self.downwind_box_endpoint = pnt
                self.downwind_box = poly
                self.downwind_box_azimuth_adjusted = azimuth
                if (not(logger is None)): logger.info("Accepted: rotation: {} average: {} azimuth: {}".format(r, average , azimuth))
                self.downwind_box_count = c
                self.shape_enhancement = s
            else:
                if (not(logger is None)): logger.info("Rejected: rotation: {} average: {} azimuth: {}".format(r, average , azimuth))
            
            self.rotations.append(Rotation(r, azimuth, average))

        if (not(logger is None)): logger.info("Forward point 0.4 lon:{}, lat:{}".format(self.downwind_box_endpoint.x, self.downwind_box_endpoint.y))
        if (not(logger is None)): logger.info("Box UL corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[0][0], self.downwind_box.exterior.coords[0][1]))
        if (not(logger is None)): logger.info("Box BL corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[1][0], self.downwind_box.exterior.coords[1][1]))
        if (not(logger is None)): logger.info("Box BR corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[2][0], self.downwind_box.exterior.coords[2][1]))
        if (not(logger is None)): logger.info("Box UR corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[3][0], self.downwind_box.exterior.coords[3][1]))

        with open(self.picklename, mode="wb") as f:
            dump(self, f)
            f.close()

        return

    def _calculate_downwindbox_enhancement_statistics(self, poly, scans, pixels, lon_bounds, lat_bounds, lons, lats, scan) -> tuple[float, float, float]:
        '''
        Local method calculating enhancement statistics.
        Note that used geometry and masking are specified by configuration values
        geometry - self.config.Algorithm_Sadavarte2021_downwindbox_geometry
        mask - self.config.Algorithm_Sadavarte2021_downwindbox_mask

        Parameters
        ----------
        poly: Polygon
            Shape within which enhancement calculation takes place.
        scans: int
            Length of the first dimension of lon_bounds, lat_bounds, lons, lats and scan arrays
        pixels: int
            Length of the second dimension of lon_bounds, lat_bounds, lons, lats and scan arrays
        lon_bounds: array[0:scans, 0:pixels, 4]
            Longitudes of TROPOMI pixels corners
        lat_bounds: array[0:scans, 0:pixels, 4]
            Latitudes of TROPOMI pixels corners
        lons: array[[0:scans, 0:pixels]
            Longitudes TROPOMI pixels centers
        lats: array[[0:scans, 0:pixels]
            Latitudes TROPOMI pixels centers
        scan: array[[0:scans, 0:pixels]

        Returns
        -------
            tuple[count, sum, average]:
                count:int - Count of pixels within polygon satisfying geometry and mask conditions. This is used by direction adjustment.
                sum:float - Sum of enhancement values. This is used by lenght and width adjustment.
                average:float - count/sum. This is used by direction adjustment.
        '''
        count = 0
        sum = 0 

        for s in range(scans):
            for p in range(pixels):
                point = Point(lons[s,p], lats[s,p])
                ul = (lon_bounds[s,p,3], lat_bounds[s,p,3])
                bl = (lon_bounds[s,p,0], lat_bounds[s,p,0])
                br = (lon_bounds[s,p,1], lat_bounds[s,p,1])
                ur = (lon_bounds[s,p,2], lat_bounds[s,p,2])
                pixelbox = Polygon([ul,bl,br,ur])
                if self.config.Algorithm_CSF_downwindbox_geometry == Geometry.CENTER and not poly.contains(point): continue
                if self.config.Algorithm_CSF_downwindbox_geometry == Geometry.CONTAINS and not poly.contains(pixelbox): continue
                if self.config.Algorithm_CSF_downwindbox_geometry == Geometry.INTERSECTS and not poly.intersects(pixelbox): continue
                
                v = scan[s,p] - self.background_enhancement 
                if self.config.Algorithm_CSF_downwindbox_mask == Mask.NEGATIVE and v < 0: continue
                sum += v
                count += 1
        return count, sum, sum/count

    def _calculate_downwindbox_lengthadjustment(self):
        '''
        From Sadavarte 2021: After we set the new wind direction, the length of the rectangular mask
        in the downwind direction (along the x-axis) is varied to define the end of the plume. This end is fixed by incrementing the
        length of the rectangular mask by 0.1° intervals until the difference between methane enhancement of two consecutive
        increments is less than 5 ppb.

        On return this method sets
            self.enhancement
            self.downwind_box_adjusted
            self.downwind_box_length
            self.downwind_endpoint
        '''
        if not (self.status == Status.PROCESSING.value):
            return
        
        if (None == self.downwind_box):
            self._calculate_downwindbox_directionadjustment()
        
        # initial enhancement was set in direction adjustment
        
        if (not(logger is None)): logger.info("Initial length: {} (m) enhancement: {} (ppb)".format(self.downwind_box_length, self.shape_enhancement))
        self.enhancement_length.append(EnhancementLength(self.downwind_box_length, self.shape_enhancement, None))
        count = 1
        newlength = self.downwind_box_length
        while True:
            newlength += Algorithm_CSF.degree_tenth
            p1 = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=self.downwind_box_azimuth_adjusted, distances=newlength)
            pnt = Point(p1[0,0], p1[0,1])
            pul = Algorithm_CSF.g.direct(points=pnt.coords, azimuths=self.downwind_box_azimuth_adjusted - 90, distances=Algorithm_CSF.degree_tenth)
            pbl = Algorithm_CSF.g.direct(points=pnt.coords, azimuths=self.downwind_box_azimuth_adjusted + 90, distances=Algorithm_CSF.degree_tenth)
            pur = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=self.downwind_box_azimuth_adjusted-90, distances=Algorithm_CSF.degree_tenth)
            pbr = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=self.downwind_box_azimuth_adjusted+90, distances=Algorithm_CSF.degree_tenth)
            poly = Polygon([(pul[0,0], pul[0,1]), (pbl[0,0], pbl[0,1]), (pbr[0,0], pbr[0,1]), (pur[0,0], pur[0,1])])
            [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, lons, lats, scan]  = self.tropomi.narrow_to_domain(self.config.Algorithm_CSF_downwindbox_geometry, poly)
            scans = maxscan - minscan
            pixels = maxpixel - minpixel
            if (scan is None): 
                if (not(logger is None)): logger.info("Rejected: length adjustment: {} Scan does not contain pixels".format(newlength))
                continue

            [c, s, _] = self._calculate_downwindbox_enhancement_statistics(poly, scans, pixels, lonbox, latbox, lons, lats, scan)

            increment = s - self.shape_enhancement
            self.enhancement_length.append(EnhancementLength(newlength, s, increment))

            if self.config.Algorithm_CSF_downwindbox_enhancement_delta < increment:
                self.shape_enhancement = s
                self.downwind_box = poly
                self.downwind_box_length = newlength
                self.downwind_box_endpoint = pnt
                if (not(logger is None)): logger.info("Accepted: Iteration: {} Length: {} (m) enhancement: {} ppb Increment: {} ".format(count, self.downwind_box_length, self.shape_enhancement, increment))
                self.downwind_box_count = c
                count += 1
            else:
                if (not(logger is None)): logger.info("Rejected: Iteration: {} Length: {} (m) enhancement: {} ppb Increment: {} Terminating length adjustment".format(count, newlength, s, increment))
                break
            

        if (not(logger is None)): logger.info("Forward point lon:{}, lat:{}".format(self.downwind_box_endpoint.x, self.downwind_box_endpoint.y))
        if (not(logger is None)): logger.info("Box BL corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[0][0], self.downwind_box.exterior.coords[0][1]))
        if (not(logger is None)): logger.info("Box UL corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[1][0], self.downwind_box.exterior.coords[1][1]))
        if (not(logger is None)): logger.info("Box UR corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[2][0], self.downwind_box.exterior.coords[2][1]))
        if (not(logger is None)): logger.info("Box BR corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[3][0], self.downwind_box.exterior.coords[3][1]))
        
        
        with open(self.picklename, mode="wb") as f:
            dump(self, f)
            f.close()
        return

    def _calculate_downwindbox_widthadjustment(self):
        '''
        From Sadavarte 2021: Similarly, the width of the rectangular mask (along the y-axis) was fixed by incrementing
        the width in the lateral direction of the plume at an interval of 0.05° until the incremental change in methane enhancement is
        less than 5 ppb.
        '''
        if not (self.status == Status.PROCESSING.value):
            return

        # This calculation depends on length adjustment to be run first
        if (None == self.shape_enhancement):
            self._calculate_downwindbox_lengthadjustment()

        if (not(logger is None)): logger.info("Initial half width: {} (m) enhancement: {} (ppb)".format(Algorithm_CSF.degree_tenth, self.shape_enhancement))
        self.enhancement_width.append(EnhancementWidth(Algorithm_CSF.degree_tenth, self.shape_enhancement, None))
        count = 1
        newwidth = Algorithm_CSF.degree_tenth
        while True:
            newwidth += Algorithm_CSF.degree_twentieth
            pul = Algorithm_CSF.g.direct(points=self.downwind_box_endpoint.coords, azimuths=self.downwind_box_azimuth_adjusted - 90, distances=newwidth)
            pbl = Algorithm_CSF.g.direct(points=self.downwind_box_endpoint.coords, azimuths=self.downwind_box_azimuth_adjusted + 90, distances=newwidth)
            pur = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=self.downwind_box_azimuth_adjusted-90, distances=newwidth)
            pbr = Algorithm_CSF.g.direct(points=self.source.xy.coords, azimuths=self.downwind_box_azimuth_adjusted+90, distances=newwidth)
            poly = Polygon([(pul[0,0], pul[0,1]), (pbl[0,0], pbl[0,1]), (pbr[0,0], pbr[0,1]), (pur[0,0], pur[0,1])])
            [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, lons, lats, scan]  = self.tropomi.narrow_to_domain(self.config.Algorithm_CSF_downwindbox_geometry, poly)
            scans = maxscan - minscan
            pixels = maxpixel - minpixel
            [c, s, _] = self._calculate_downwindbox_enhancement_statistics(poly, scans, pixels, lonbox, latbox, lons, lats, scan)

            increment = s - self.shape_enhancement
            self.enhancement_width.append(EnhancementWidth(newwidth, s, increment))

            if self.config.Algorithm_CSF_downwindbox_enhancement_delta < increment:
                self.shape_enhancement = s
                self.downwind_box = poly
                self.downwind_box_half_width = newwidth
                if (not(logger is None)): logger.info("Accepted Iteration: {} half-width: {} (m) enhancement: {} ppb Increment: {} ".format(count, self.downwind_box_half_width, self.shape_enhancement, increment))
                self.downwind_box_count = c
                count += 1
            else:
                if (not(logger is None)): logger.info("Rejected Iteration: {} half-width: {} (m) enhancement: {} ppb Increment: {}. Terminating width adjustment ".format(count, newwidth, s, increment))
                break

        if (not(logger is None)): logger.info("Box BR corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[0][0], self.downwind_box.exterior.coords[0][1]))
        if (not(logger is None)): logger.info("Box TR corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[1][0], self.downwind_box.exterior.coords[1][1]))
        if (not(logger is None)): logger.info("Box TL corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[2][0], self.downwind_box.exterior.coords[2][1]))
        if (not(logger is None)): logger.info("Box BL corner lon:{}, lat:{}".format(self.downwind_box.exterior.coords[3][0], self.downwind_box.exterior.coords[3][1]))
        
        
        with open(self.picklename, mode="wb") as f:
            dump(self, f)
            f.close()

        return
  
    def _chart_downwindbox_imagecoords_positive(self, fileout:str):
        '''
        Chart positive pixels in downwind box using image coordinates (scan,pixel)
        '''
        if (0 < self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value) or (0 < self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            print("Failure: Algorithm failed with status: {}".format(self.status))
            return
        
        title = "Algorithm CSF. Positive. Orbit: {} Processor: {}".format(self.tropomi.orbit, self.tropomi.processor_version)

        [minscan, minpixel, maxscan, maxpixel, _, _, _, _, _] = self.tropomi.narrow_to_domain(Geometry.INTERSECTS, self.downwind_box)
        scanCH4 = self.tropomi.get_variable_data("methane_mixing_ratio_bias_corrected")

        scan1 = scanCH4 - self.background

        for s in range(minscan, maxscan):
            for p in range(minpixel, maxpixel):
                if scan1[0,s,p] < 0: scan1[0,s,p] = ma.masked

        count = ma.count(scan1)
        frac_y_x = (maxpixel - minpixel) / (maxscan - minscan)
        frac_y_x = 10 * frac_y_x + 0.2 # 0.2 correct for text under

        v_min = 0
        v_max = 0
        
        if 0 < count:
            v_min = ma.min(scan1[0, minscan:maxscan, minpixel:maxpixel])
            v_max = ma.max(scan1[0, minscan:maxscan, minpixel:maxpixel])
        
        if (not(logger is None)): 
            logger.info("Scan count: {} min: {} max: {}".format(count, v_min, v_max))

        fig = plt.figure(figsize=(10,frac_y_x))
        ax = plt.axes()
        ax.set_xlim(left=minpixel,  right=maxpixel + 1)
        ax.set_ylim(bottom=minscan,  top=maxscan + 1)
        
        # Draw grid of TROPOMI data
        if 0 < count:
            cs = ax.pcolormesh(scan1[0,:,:], shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.rainbow)

            # Draw value bar below chart
            cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
            cbar.ax.set_xlabel("CH4 ppb")

        # Draw mine marker in image coordinates
        ax.plot(self.tropomi_source_pixel, self.tropomi_source_scan, 'go', markersize=7)
        ax.text(self.tropomi_source_pixel + 0.25, self.tropomi_source_scan + 0.25, self.source.display_name, color="black")

        plt.title(title)

        if self.config.Algorithm_CSF_transect_background == "algorithm":
            if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
                text3 = "Background: Geometry: {}, Algorithm: Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
            else:
                text3 = "Background: Geometry: {}, Algorithm: Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "domain":
            text3 = "Background: Geometry: {}, Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "manual" and not (self.config.Algorithm_CSF_transect_background_value is None):
            text3 = "Background: Manually set {:.3f} (ppb) ".format(self.background)
        elif self.config.Algorithm_CSF_transect_background == "minimum":
            text3 = "Background: Geometry: {}, Minimum {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "upwind":
            text3 = "Background: Geometry: {}, Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)

        text2 = "Pressure averaged boundary layer wind: u: {:.3f} (m/s) v:{:.3f} (m/s)".format(self.u, self.v)

        if self.status == Status.SUCCESS.value:
            text1 = "Sucess Q {:.3f} (t/hr)".format(self.q_valid_hour)
        elif 0 < (self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            text1 = "Fail: No data in domain"
        elif 0 < (self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value):
            text1 = "Fail: No data in downwind box"
        elif 0 < (self.status & Status.FAILURE_TRANSECTS_BELOW_3.value):
            text1 = "Fail: Number of valid transects < 3"
        elif 0 < (self.status & Status.FAILURE_WIND_BELOW2.value):
            text1 = "Fail: Wind under 2 (m/s)"

        text4 = "Enhancement: Geometry: {}, Masking: {}".format(self.config.Algorithm_CSF_downwindbox_geometry.name, self.config.Algorithm_CSF_downwindbox_mask.name)
        text5 = "Transect: Valid count: {} Minimum: {:.3f} (ppb)".format(self.transects_valid_count, self.downwind_box_minimumvalue)
        
        fig.text(0.1, 0.15, text1, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.125, text2, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.1, text3, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.075, text4, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.05, text5, horizontalalignment="left", wrap=False) 
        fig.subplots_adjust(bottom=0.2)

        if not fileout is None:
            plt.savefig(fileout)
            print("Chart generated: {}".format(fileout))
        else:
            plt.show()
        plt.close("all")

    def _chart_downwindbox_imagecoords_valid(self, filename: str):
        '''
        Chart valid pixels in downwind box using image coordinates (scan,pixel)
        '''
        if (0 < self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value) or (0 < self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            print("Failure: Algorithm failed with status: {}".format(self.status))
            return
        
        title = "Algorithm CSF. Valid. Orbit: {} Processor: {}".format(self.tropomi.orbit, self.tropomi.processor_version)

        # Takes the largest box
        [minscan, minpixel, maxscan, maxpixel, _, _, _, _, _] = self.tropomi.narrow_to_domain(Geometry.INTERSECTS, self.downwind_box)
        scanCH4 = self.tropomi.get_variable_data("methane_mixing_ratio_bias_corrected") 
        scan1 = scanCH4 - self.background

        count = ma.count(scan1)
        v_min = 0
        v_max = 0
        frac_y_x = (maxpixel - minpixel) / (maxscan - minscan)
        frac_y_x = 10 * frac_y_x + 0.2 # 0.2 correct for text under

        if 0 < count:
            v_min = ma.min(scan1[0, minscan:maxscan, minpixel:maxpixel])
            v_max = ma.max(scan1[0, minscan:maxscan, minpixel:maxpixel])
        
        if (not(logger is None)): 
            logger.info("Scan count: {} min: {} max: {}".format(count, v_min, v_max))

        fig = plt.figure(figsize=(10,frac_y_x))
        ax = plt.axes()
        ax.set_xlim(left=minpixel,  right=maxpixel + 1)
        ax.set_ylim(bottom=minscan,  top=maxscan + 1)
        
        # Draw grid of TROPOMI data
        if 0 < count:
            cs = ax.pcolormesh(scan1[0,:,:], shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.rainbow)

            # Draw value bar below chart
            cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
            cbar.ax.set_xlabel("CH4 ppb")

        # Draw mine marker in image coordinates
        ax.plot(self.tropomi_source_pixel, self.tropomi_source_scan, 'go', markersize=7)
        ax.text(self.tropomi_source_pixel + 0.25, self.tropomi_source_scan + 0.25, self.source.display_name, color="black")

        plt.title(title)

        if self.config.Algorithm_CSF_transect_background == "algorithm":
            if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
                text3 = "Background: Geometry: {}, Algorithm: Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
            else:
                text3 = "Background: Geometry: {}, Algorithm: Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "domain":
            text3 = "Background: Geometry: {}, Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "manual" and not (self.config.Algorithm_CSF_transect_background_value is None):
            text3 = "Background: Manually set {:.3f} (ppb) ".format(self.background)
        elif self.config.Algorithm_CSF_transect_background == "minimum":
            text3 = "Background: Geometry: {}, Minimum {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "upwind":
            text3 = "Background: Geometry: {}, Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)

        text2 = "Pressure averaged boundary layer wind: u: {:.3f} (m/s) v:{:.3f} (m/s)".format(self.u, self.v)

        if self.status == Status.SUCCESS.value:
            text1 = "Sucess Q {:.3f} (t/hr)".format(self.q_valid_hour)
        elif 0 < (self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            text1 = "Fail: No data in domain"
        elif 0 < (self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value):
            text1 = "Fail: No data in downwind box"
        elif 0 < (self.status & Status.FAILURE_TRANSECTS_BELOW_3.value):
            text1 = "Fail: Number of valid transects < 3"
        elif 0 < (self.status & Status.FAILURE_WIND_BELOW2.value):
            text1 = "Fail: Wind under 2 (m/s)"

        text4 = "Enhancement: Geometry: {}, Masking: {}".format(self.config.Algorithm_CSF_downwindbox_geometry.name, self.config.Algorithm_CSF_downwindbox_mask.name)
        text5 = "Transect: Valid count: {} Minimum: {:.3f} (ppb)".format(self.transects_valid_count, self.downwind_box_minimumvalue)
        
        fig.text(0.1, 0.15, text1, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.125, text2, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.1, text3, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.075, text4, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.05, text5, horizontalalignment="left", wrap=False) 
        fig.subplots_adjust(bottom=0.2)

        if not filename is None:
            plt.savefig(filename)
            print("Chart generated {}".format(filename))
        else:
            plt.show()
        plt.close("all")
        pass

    def _chart_combined(self, filename: str):
        '''
        Chart valid and positive graphs on the same chart
        '''
        # Data    
        multipolygon = self.downwind_box
        iscombined = False
        if ((self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count)): 
            multipolygon = self.downwind_box.union(self.upwind_box)
            iscombined = True

        (minlon, minlat, maxlon, maxlat) = multipolygon.bounds

        # Takes the largest box
        [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, lons, lats, scanCH4] = self.tropomi.narrow_to_domain(Geometry.INTERSECTS, multipolygon)
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        
        scan1 = scanCH4 - self.background
        scan2 = scanCH4 - self.background

        for s in range(scans):
            for p in range(pixels):
                point = Point(lons[s,p], lats[s,p])
                ul = (lonbox[s,p,3], latbox[s,p,3])
                bl = (lonbox[s,p,0], latbox[s,p,0])
                br = (lonbox[s,p,1], latbox[s,p,1])
                ur = (lonbox[s,p,2], latbox[s,p,2])
                pixelbox = Polygon([ul,bl,br,ur])

                # We are dealing with multi polygon so deal with pixels intersecting the upwind box
                if iscombined and self.upwind_box.intersects(pixelbox):
                    if (self.config.Algorithm_CSF_background_geometry == Geometry.CENTER) and not self.upwind_box.contains(point): 
                        scan1[s,p] = ma.masked
                        scan2[s,p] = ma.masked
                    if (self.config.Algorithm_CSF_background_geometry == Geometry.CONTAINS) and not self.upwind_box.contains(pixelbox): 
                        scan1[s,p] = ma.masked
                        scan2[s,p] = ma.masked
                    
                # We are dealing with multi polygon so deal with pixels intersecting the downwind box
                elif self.downwind_box.intersects(pixelbox): 
                    if scan2[s,p] < 0: scan2[s,p] = ma.masked
                
                # Pixels outside upwind_box and donwind_box 
                else: 
                    scan1[s,p] = ma.masked
                    scan2[s,p] = ma.masked

        count = ma.count(scan1)
        v_min = 0
        v_max = 0
        if 0 < count:
            v_min = ma.min(scan1)
            v_max = ma.max(scan1)

        if (not(logger is None)): 
            logger.info("Scan count: {} min: {} max: {}".format(count, v_min, v_max))

        # Charting
        plateCarree = ccrs.PlateCarree()
        geodetic = ccrs.Geodetic()
        fig = plt.figure(figsize=(10,5), layout=None)
        fig.suptitle("Algorithm CSF. Valid. Orbit: {} Processor: {}".format(self.tropomi.orbit, self.tropomi.processor_version))
        gs = GridSpec(3,2, left=0.1, right=0.9)
        axs00 = fig.add_subplot(gs[:-1, 0], projection=plateCarree)
        axs01 = fig.add_subplot(gs[:-1, 1], projection=plateCarree)
        axs10 = fig.add_subplot(gs[2,:])
        axs10.set_axis_off()

        axs00.set_extent(extents=[minlon-0.1, maxlon+0.1, minlat-0.1, maxlat+0.1], crs=plateCarree)
        axs00.coastlines()

        axs01.set_extent(extents=[minlon-0.1, maxlon+0.1, minlat-0.1, maxlat+0.1], crs=plateCarree)
        axs01.coastlines()
        
        # Draw grid of TROPOMI data
        if 0 < count:
            cs1 = axs00.pcolormesh(lons, lats, scan1, shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.rainbow, transform=plateCarree)
            cs2 = axs01.pcolormesh(lons, lats, scan2, shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.rainbow, transform=plateCarree)
            cbar1 = fig.colorbar(cs1, ax=[axs00, axs01] , orientation="horizontal", location="bottom", pad=0.1, shrink=0.6)
            cbar1.set_label("CH4 ppb")

        # Draw wind arrow at location of ERA 5 retrieval
        axs00.quiver(self.source.gridpoint.x, self.source.gridpoint.y, self.u, self.v)
        axs01.quiver(self.source.gridpoint.x, self.source.gridpoint.y, self.u, self.v)

        # Draw mine marker
        Source.plot_source_in_extent(self.tropomi_source_date,minlon, minlat, maxlon, maxlat, axs00)
        Source.plot_source_in_extent(self.tropomi_source_date,minlon, minlat, maxlon, maxlat, axs01)

        # Plot point at which wind is calculated
        axs00.plot(self.source.gridpoint.x, self.source.gridpoint.y, 'gv', markersize=7)
        axs01.plot(self.source.gridpoint.x, self.source.gridpoint.y, 'gv', markersize=7)

        # Draw upwind box if it is used
        if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
            drawTrajectory(axs00, geodetic, plateCarree, self.upwind_box.exterior.coords)        
            drawTrajectory(axs01, geodetic, plateCarree, self.upwind_box.exterior.coords)        

        # Draw plume shape
        drawTrajectory(axs00, geodetic, plateCarree, self.downwind_box.exterior.coords)
        drawTrajectory(axs01, geodetic, plateCarree, self.downwind_box.exterior.coords)

        # Draw transects
        for t in self.transects:
            drawTrajectory(axs00, geodetic, plateCarree, t.line.coords)
            drawTrajectory(axs01, geodetic, plateCarree, t.line.coords)

        gl1 = axs00.gridlines(crs=plateCarree, draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl1.top_labels = False
        gl1.left_labels = False
        gl1.right_labels = False
        gl1.xformatter = LONGITUDE_FORMATTER
        gl1.yformatter = LATITUDE_FORMATTER
        gl1.xlabel_style = {'size': 10, 'color': 'gray'}
        gl1.ylabel_style = {'size': 10, 'color': 'gray'}

        gl2 = axs01.gridlines(crs=plateCarree, draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl2.top_labels = False
        gl2.left_labels = False
        gl2.xformatter = LONGITUDE_FORMATTER
        gl2.yformatter = LATITUDE_FORMATTER
        gl2.xlabel_style = {'size': 10, 'color': 'gray'}
        gl2.ylabel_style = {'size': 10, 'color': 'gray'}

        axs00.set_title("Valid")
        axs01.set_title("Positive")
        if self.config.Algorithm_CSF_transect_background == "algorithm":
            if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
                text3 = "Background: Geometry: {}, Algorithm: Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
            else:
                text3 = "Background: Geometry: {}, Algorithm: Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "domain":
            text3 = "Background: Geometry: {}, Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "manual" and not (self.config.Algorithm_CSF_transect_background_value is None):
            text3 = "Background: Manually set {:.3f} (ppb) ".format(self.background)
        elif self.config.Algorithm_CSF_transect_background == "minimum":
            text3 = "Background: Geometry: {}, Minimum {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "upwind":
            text3 = "Background: Geometry: {}, Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)

        text2 = "Pressure averaged boundary layer wind: u: {:.3f} (m/s) v:{:.3f} (m/s)".format(self.u, self.v)

        if self.status == Status.SUCCESS.value:
            text1 = "Sucess Q {:.3f} (t/hr)".format(self.q_valid_hour)
        elif 0 < (self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            text1 = "Fail: No data in domain"
        elif 0 < (self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value):
            text1 = "Fail: No data in downwind box"
        elif 0 < (self.status & Status.FAILURE_TRANSECTS_BELOW_3.value):
            text1 = "Fail: Number of valid transects < 3"
        elif 0 < (self.status & Status.FAILURE_WIND_BELOW2.value):
            text1 = "Fail: Wind under 2 (m/s)"

        text4 = "Enhancement: Geometry: {}, Masking: {}".format(self.config.Algorithm_CSF_downwindbox_geometry.name, self.config.Algorithm_CSF_downwindbox_mask.name)
        text5 = "Transect: Valid count: {} Minimum: {:.3f} (ppb)".format(self.transects_valid_count, self.downwind_box_minimumvalue)
        
        axs10.text(0, 0.6, text1, horizontalalignment="left", wrap=False) 
        axs10.text(0, 0.45, text2, horizontalalignment="left", wrap=False) 
        axs10.text(0, 0.3, text3, horizontalalignment="left", wrap=False) 
        axs10.text(0, 0.15, text4, horizontalalignment="left", wrap=False) 
        axs10.text(0, 0, text5, horizontalalignment="left", wrap=False) 

        if not filename is None:
            plt.savefig(filename)
            print("Chart generated {}".format(filename)) 
        else:
            plt.show()
        plt.close("all")

    def chart_lonlat_all_elements(self, filename: str):
        '''
        Chart all the initial downwind box, valid pixels, and upwind box
        '''
        if (0 < self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value) or (0 < self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            print("Failure: Algorithm failed with status: {}".format(self.status))
            return
        
        title = "Algorithm CSF. Valid. Orbit: {} Processor: {}".format(self.tropomi.orbit, self.tropomi.processor_version)

        multipolygon = self.downwind_box
        iscombined = False
        if ((self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count)): 
            multipolygon = self.downwind_box.union(self.upwind_box)
            iscombined = True

        (minlon, minlat, maxlon, maxlat) = multipolygon.bounds

        # Takes the largest box
        [minscan, minpixel, maxscan, maxpixel,  lonbox, latbox, lons, lats, scanCH4] = self.tropomi.narrow_to_domain(Geometry.INTERSECTS, multipolygon)
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        scan1 = scanCH4 - self.background

        for s in range(scans):
            for p in range(pixels):
                point = Point(lons[s,p], lats[s,p])
                ul = (lonbox[s,p,3], latbox[s,p,3])
                bl = (lonbox[s,p,0], latbox[s,p,0])
                br = (lonbox[s,p,1], latbox[s,p,1])
                ur = (lonbox[s,p,2], latbox[s,p,2])
                pixelbox = Polygon([ul,bl,br,ur])

                # We are dealing with multi polygon so deal with pixels intersecting the upwind box
                if iscombined and self.upwind_box.intersects(pixelbox):
                    if (self.config.Algorithm_CSF_background_geometry == Geometry.CENTER) and not self.upwind_box.contains(point): scan1[s,p] = ma.masked
                    if (self.config.Algorithm_CSF_background_geometry == Geometry.CONTAINS) and not self.upwind_box.contains(pixelbox): scan1[s,p] = ma.masked
                    
                # We are dealing with multi polygon so deal with pixels intersecting the downwind box
                elif self.downwind_box.intersects(pixelbox): continue
                # Pixels outside upwind_box and donwind_box 
                else: scan1[s,p] = ma.masked

        count = ma.count(scan1)
        minlon = ma.min(lons)
        maxlon = ma.max(lons)
        minlat = ma.min(lats)
        maxlat = ma.max(lats)

        frac_y_x = (maxlat - minlat) / (maxlon - minlon)
        frac_y_x = 10 * frac_y_x + 0.2 # 0.2 correct for text under

        v_min = 0
        v_max = 0
        if 0 < count:
            v_min = ma.min(scan1)
            v_max = ma.max(scan1)
        
        if (not(logger is None)): 
            logger.info("Scan count: {} min: {} max: {}".format(count, v_min, v_max))

        plateCarree = ccrs.PlateCarree()
        geodetic = ccrs.Geodetic()
        fig = plt.figure(figsize=(10,frac_y_x))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extents=[minlon-0.1, maxlon+0.1, minlat-0.1, maxlat+0.1], crs=plateCarree)
        ax.coastlines()
        
        # Draw grid of TROPOMI data
        if 0 < count:
            cs = ax.pcolormesh(lons, lats, scan1, shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.rainbow, transform=ccrs.PlateCarree())

            # Draw value bar below chart
            cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
            cbar.ax.set_xlabel("CH4 ppb")

        # Draw wind arrow at location of ERA 5 retrieval
        ax.quiver(self.source.gridpoint.x, self.source.gridpoint.y, self.u, self.v)


        # Draw mine marker
        Source.plot_source_in_extent(self.tropomi_source_date, minlon, minlat, maxlon, maxlat, ax)

        # Plot point at which wind is calculated
        ax.plot(self.source.gridpoint.x, self.source.gridpoint.y, 'gv', markersize=7)

        # Draw upwind box if it is used
        if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
            drawTrajectory(plt, geodetic, plateCarree, self.upwind_box.exterior.coords)        

        # Draw the intial downwind box
        drawTrajectory(plt, geodetic, plateCarree, self.downwind_box_initial.exterior.coords)

        # Draw plume shape
        drawTrajectory(plt, geodetic, plateCarree, self.downwind_box.exterior.coords)

        # Draw transects
        for t in self.transects:
            drawTrajectory(plt, geodetic, plateCarree, t.line.coords)

        plt.title(title)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'}

        if self.config.Algorithm_CSF_transect_background == "algorithm":
            if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
                text3 = "Background: Geometry: {}, Algorithm: Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
            else:
                text3 = "Background: Geometry: {}, Algorithm: Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "domain":
            text3 = "Background: Geometry: {}, Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "manual" and not (self.config.Algorithm_CSF_transect_background_value is None):
            text3 = "Background: Manually set {:.3f} (ppb) ".format(self.background)
        elif self.config.Algorithm_CSF_transect_background == "minimum":
            text3 = "Background: Geometry: {}, Minimum {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "upwind":
            text3 = "Background: Geometry: {}, Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)

        text2 = "Pressure averaged boundary layer wind: u: {:.3f} (m/s) v:{:.3f} (m/s)".format(self.u, self.v)

        if self.status == Status.SUCCESS.value:
            text1 = "Sucess Q {:.3f} (t/hr)".format(self.q_valid_hour)
        elif 0 < (self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            text1 = "Fail: No data in domain"
        elif 0 < (self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value):
            text1 = "Fail: No data in downwind box"
        elif 0 < (self.status & Status.FAILURE_TRANSECTS_BELOW_3.value):
            text1 = "Fail: Number of valid transects < 3"
        elif 0 < (self.status & Status.FAILURE_WIND_BELOW2.value):
            text1 = "Fail: Wind under 2 (m/s)"

        text4 = "Enhancement: Geometry: {}, Masking: {}".format(self.config.Algorithm_CSF_downwindbox_geometry.name, self.config.Algorithm_CSF_downwindbox_mask.name)
        text5 = "Transect: Valid count: {} Minimum: {:.3f} (ppb)".format(self.transects_valid_count, self.downwind_box_minimumvalue)
        
        fig.text(0.1, 0.15, text1, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.125, text2, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.1, text3, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.075, text4, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.05, text5, horizontalalignment="left", wrap=False) 
        fig.subplots_adjust(bottom=0.2)

        if not filename is None:
            plt.savefig(filename)            
            print("Chart generated: {}".format(filename))
        else:
            plt.show()
        plt.close("all")

    def chart_lonlat_valid(self, filename: str):
        '''
        Chart all valid pixels replicating figure S1 in Sadavarte 2021
        '''
        if (0 < self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value) or (0 < self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            print("Failure: Algorithm failed with status: {}".format(self.status))
            return
        
        title = "Algorithm CSF. Valid. Orbit: {} Processor: {}".format(self.tropomi.orbit, self.tropomi.processor_version)

        multipolygon = self.downwind_box
        iscombined = False
        if (not self.upwind_box_background_count is None and (self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count)): 
            multipolygon = self.downwind_box.union(self.upwind_box)
            iscombined = True

        (minlon, minlat, maxlon, maxlat) = multipolygon.bounds

        # Takes the largest box
        [minscan, minpixel, maxscan, maxpixel,  lonbox, latbox, lons, lats, scanCH4] = self.tropomi.narrow_to_domain(Geometry.INTERSECTS, multipolygon)
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        scan1 = scanCH4 - self.background

        for s in range(scans):
            for p in range(pixels):
                point = Point(lons[s,p], lats[s,p])
                ul = (lonbox[s,p,3], latbox[s,p,3])
                bl = (lonbox[s,p,0], latbox[s,p,0])
                br = (lonbox[s,p,1], latbox[s,p,1])
                ur = (lonbox[s,p,2], latbox[s,p,2])
                pixelbox = Polygon([ul,bl,br,ur])

                # We are dealing with multi polygon so deal with pixels intersecting the upwind box
                if iscombined and self.upwind_box.intersects(pixelbox):
                    if (self.config.Algorithm_CSF_background_geometry == Geometry.CENTER) and not self.upwind_box.contains(point): scan1[s,p] = ma.masked
                    if (self.config.Algorithm_CSF_background_geometry == Geometry.CONTAINS) and not self.upwind_box.contains(pixelbox): scan1[s,p] = ma.masked
                    
                # We are dealing with multi polygon so deal with pixels intersecting the downwind box
                elif self.downwind_box.intersects(pixelbox): continue
                # Pixels outside upwind_box and donwind_box 
                else: scan1[s,p] = ma.masked

        count = ma.count(scan1)
        minlon = ma.min(lons)
        maxlon = ma.max(lons)
        minlat = ma.min(lats)
        maxlat = ma.max(lats)

        frac_y_x = (maxlat - minlat) / (maxlon - minlon)
        frac_y_x = 10 * frac_y_x + 0.2 # 0.2 correct for text under

        v_min = 0
        v_max = 0
        if 0 < count:
            v_min = ma.min(scan1)
            v_max = ma.max(scan1)
        
        if (not(logger is None)): 
            logger.info("Scan count: {} min: {} max: {}".format(count, v_min, v_max))

        plateCarree = ccrs.PlateCarree()
        geodetic = ccrs.Geodetic()
        fig = plt.figure(figsize=(10,frac_y_x))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extents=[minlon-0.1, maxlon+0.1, minlat-0.1, maxlat+0.1], crs=plateCarree)
        ax.coastlines()
        
        # Draw grid of TROPOMI data
        if 0 < count:
            cs = ax.pcolormesh(lons, lats, scan1, shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.rainbow, transform=ccrs.PlateCarree())

            # Draw value bar below chart
            cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
            cbar.ax.set_xlabel("CH4 ppb")

        # Draw wind arrow at location of ERA 5 retrieval
        ax.quiver(self.source.gridpoint.x, self.source.gridpoint.y, self.u, self.v)


        # Draw mine marker
        Source.plot_source_in_extent(self.tropomi_source_date, minlon, minlat, maxlon, maxlat, ax)

        # Plot point at which wind is calculated
        ax.plot(self.source.gridpoint.x, self.source.gridpoint.y, 'gv', markersize=7)

        # Draw upwind box if it is used
        if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
            drawTrajectory(plt, geodetic, plateCarree, self.upwind_box.exterior.coords)        

        # Draw plume shape
        drawTrajectory(plt, geodetic, plateCarree, self.downwind_box.exterior.coords)

        # Draw transects
        for t in self.transects:
            drawTrajectory(plt, geodetic, plateCarree, t.line.coords)

        plt.title(title)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'}

        if self.config.Algorithm_CSF_transect_background == "algorithm":
            if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
                text3 = "Background: Geometry: {}, Algorithm: Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
            else:
                text3 = "Background: Geometry: {}, Algorithm: Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "domain":
            text3 = "Background: Geometry: {}, Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "manual" and not (self.config.Algorithm_CSF_transect_background_value is None):
            text3 = "Background: Manually set {:.3f} (ppb) ".format(self.background)
        elif self.config.Algorithm_CSF_transect_background == "minimum":
            text3 = "Background: Geometry: {}, Minimum {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "upwind":
            text3 = "Background: Geometry: {}, Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)

        text2 = "Pressure averaged boundary layer wind: u: {:.3f} (m/s) v:{:.3f} (m/s)".format(self.u, self.v)

        if self.status == Status.SUCCESS.value:
            text1 = "Sucess Q {:.3f} (t/hr)".format(self.q_valid_hour)
        elif 0 < (self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            text1 = "Fail: No data in domain"
        elif 0 < (self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value):
            text1 = "Fail: No data in downwind box"
        elif 0 < (self.status & Status.FAILURE_TRANSECTS_BELOW_3.value):
            text1 = "Fail: Number of valid transects < 3"
        elif 0 < (self.status & Status.FAILURE_WIND_BELOW2.value):
            text1 = "Fail: Wind under 2 (m/s)"

        text4 = "Enhancement: Geometry: {}, Masking: {}".format(self.config.Algorithm_CSF_downwindbox_geometry.name, self.config.Algorithm_CSF_downwindbox_mask.name)
        text5 = "Transect: Valid count: {} Minimum: {:.3f} (ppb)".format(self.transects_valid_count, self.downwind_box_minimumvalue)
        
        fig.text(0.1, 0.15, text1, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.125, text2, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.1, text3, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.075, text4, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.05, text5, horizontalalignment="left", wrap=False) 
        fig.subplots_adjust(bottom=0.2)

        if not filename is None:
            plt.savefig(filename)            
            print("Chart generated: {}".format(filename))
        else:
            plt.show()
        plt.close("all")

    def chart_lonlat_positive(self, filename: str):
        '''
        Chart positive pixels only.
        '''
        if (0 < self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value) or (0 < self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            print("Failure: Algorithm failed with status: {}".format(self.status))
            return
        
        title = "Algorithm CSF. Positive. Orbit: {} Processor: {}".format(self.tropomi.orbit, self.tropomi.processor_version)

        multipolygon = self.downwind_box
        iscombined = False
        if (self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count): 
            multipolygon = self.downwind_box.union(self.upwind_box)
            iscombined = True

        (minlon, minlat, maxlon, maxlat) = multipolygon.bounds

        # Takes the largest box
        [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, lons, lats, scanCH4] = self.tropomi.narrow_to_domain(Geometry.INTERSECTS, multipolygon)
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        
        scan1 = scanCH4 - self.background

        for s in range(scans):
            for p in range(pixels):
                point = Point(lons[s,p], lats[s,p])
                ul = (lonbox[s,p,3], latbox[s,p,3])
                bl = (lonbox[s,p,0], latbox[s,p,0])
                br = (lonbox[s,p,1], latbox[s,p,1])
                ur = (lonbox[s,p,2], latbox[s,p,2])
                pixelbox = Polygon([ul,bl,br,ur])

                # We are dealing with multi polygon so deal with pixels intersecting the upwind box
                if iscombined and self.upwind_box.intersects(pixelbox):
                    if (self.config.Algorithm_CSF_background_geometry == Geometry.CENTER) and not self.upwind_box.contains(point): scan1[s,p] = ma.masked
                    if (self.config.Algorithm_CSF_background_geometry == Geometry.CONTAINS) and not self.upwind_box.contains(pixelbox): scan1[s,p] = ma.masked
                    
                # We are dealing with multi polygon so deal with pixels intersecting the downwind box
                elif self.downwind_box.intersects(pixelbox): 
                    if scan1[s,p] < 0: scan1[s,p] = ma.masked

                # Pixels outside upwind_box and donwind_box 
                else: scan1[s,p] = ma.masked

        count = ma.count(scan1)
        minlon = ma.min(lons)
        maxlon = ma.max(lons)
        minlat = ma.min(lats)
        maxlat = ma.max(lats)

        frac_y_x = (maxlat - minlat) / (maxlon - minlon)
        frac_y_x = 10 * frac_y_x + 0.2 # 0.2 correct for text under

        v_min = 0
        v_max = 0
        if 0 < count:
            v_min = ma.min(scan1)
            v_max = ma.max(scan1)
        
        if (not(logger is None)): 
            logger.info("Scan count: {} min: {} max: {}".format(count, v_min, v_max))

        plateCarree = ccrs.PlateCarree()
        geodetic = ccrs.Geodetic()
        fig = plt.figure(figsize=(10, frac_y_x))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extents=[minlon-0.1, maxlon+0.1, minlat-0.1, maxlat+0.1], crs=plateCarree)
        ax.coastlines()
        
        # Draw grid of TROPOMI data
        if 0 < count:
            cs = ax.pcolormesh(lons, lats, scan1, shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.rainbow, transform=ccrs.PlateCarree())

            # Draw value bar below chart
            # cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
            cbar = fig.colorbar(cs, orientation="horizontal")
            cbar.ax.set_xlabel("CH4 ppb")

        # Draw wind arrow at location of ERA 5 retrieval
        ax.quiver(self.source.gridpoint.x, self.source.gridpoint.y, self.u, self.v)

        # Draw mine marker
        Source.plot_source_in_extent(self.tropomi_source_date, minlon, minlat, maxlon, maxlat, ax)

        # Plot point at which wind is calculated
        ax.plot(self.source.gridpoint.x, self.source.gridpoint.y, 'gv', markersize=7)

        # Draw upwind box if it is used
        if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
            drawTrajectory(plt, geodetic, plateCarree, self.upwind_box.exterior.coords)        

        # Draw plume shape
        drawTrajectory(plt, geodetic, plateCarree, self.downwind_box.exterior.coords)

        # Draw transects
        for t in self.transects:
            drawTrajectory(plt, geodetic, plateCarree, t.line.coords)

        plt.title(title)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'}

        if self.config.Algorithm_CSF_transect_background == "algorithm":
            if self.config.Algorithm_CSF_background_minimum_count < self.upwind_box_background_count:
                text3 = "Background: Geometry: {}, Algorithm: Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
            else:
                text3 = "Background: Geometry: {}, Algorithm: Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "domain":
            text3 = "Background: Geometry: {}, Domain median {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)
        elif self.config.Algorithm_CSF_transect_background == "manual" and not (self.config.Algorithm_CSF_transect_background_value is None):
            text3 = "Background: Manually set {:.3f} (ppb) ".format(self.background)
        elif self.config.Algorithm_CSF_transect_background == "minimum":
            text3 = "Background: Minimum {:.3f} (ppb) ".format(self.background)
        elif self.config.Algorithm_CSF_transect_background == "upwind":
            text3 = "Background: Geometry: {}, Upwind average {:.3f} (ppb) ".format(self.config.Algorithm_CSF_background_geometry.name, self.background)

        text2 = "Pressure averaged boundary layer wind: u: {:.3f} (m/s) v:{:.3f} (m/s)".format(self.u, self.v)

        if self.status == Status.SUCCESS.value:
            text1 = "Sucess Q {:.3f} (t/hr)".format(self.q_valid_hour)
        elif 0 < (self.status & Status.FAILURE_NODATA_INDOMAIN.value):
            text1 = "Fail: No data in domain"
        elif 0 < (self.status & Status.FAILURE_NODATA_DOWNWINDBOX.value):
            text1 = "Fail: No data in downwind box"
        elif 0 < (self.status & Status.FAILURE_TRANSECTS_BELOW_3.value):
            text1 = "Fail: Number of valid transects < 3"
        elif 0 < (self.status & Status.FAILURE_WIND_BELOW2.value):
            text1 = "Fail: Wind under 2 (m/s)"

        text4 = "Enhancement: Geometry: {}, Masking: {}".format(self.config.Algorithm_CSF_downwindbox_geometry.name, self.config.Algorithm_CSF_downwindbox_mask.name)
        text5 = "Transect: Valid count: {} Minimum: {:.3f} (ppb)".format(self.transects_valid_count, self.downwind_box_minimumvalue)
        
        fig.text(0.1, 0.15, text1, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.125, text2, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.1, text3, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.075, text4, horizontalalignment="left", wrap=False) 
        fig.text(0.1, 0.05, text5, horizontalalignment="left", wrap=False) 
        fig.subplots_adjust(bottom=0.2)

        if not filename is None:
            plt.savefig(filename)
            print("Chart generated: {}".format(filename))
        else:
            plt.show()
        plt.close("all")

    def chart_transects(self, filename: str) -> None:
        '''
        Chart emission rates for valid transects
        '''
        windspeed = sqrt(self.u * self.u + self.v * self.v) # (m/s)
        direction_correction = cos(pi * (self.downwind_box_azimuth_adjusted - self.downwind_box_azimuth_initial) / 180.0)
        x = []
        y = []
        i = 0
        for t in self.transects:
            if t.isvalid():
                x.append(i + 4)
                y.append(t.enhkg * 3.6 * windspeed * direction_correction) # convert kg/s to t/hour
            i += 1
        plt.plot(x, y, marker="x", linestyle="")
        plt.grid(axis="y", color="gray", linestyle="--", linewidth=0.5 )
        plt.title("CSF. Transects. Orbit: {} Processor: {}.".format(self.tropomi.orbit, self.tropomi.processor_version))
        plt.xlabel("Transect")
        plt.ylabel("Emission rate (t/hour)")
        if not filename is None:
            plt.savefig(filename)
            print("Chart generated: {}".format(filename))
        else:
            plt.show()
        plt.close()
        return
    
    def calculate_sourcerate(self):
        '''
        From Sadavarte 2021: For a daily source rate, we take the mean of all of the emission estimates calculated for individual transects (j = 1,..., n,
        where n is the number of transects) between the source and the end of the plume. 
        ...
        With this requirement, we only calculate the source rate from plumes with at least three or more transects.
        '''
        if not (self.status == Status.PROCESSING.value):
            # move here handlers for terminal errors which stopped calculations on the way
            if (not(logger is None)): logger.info("FAILURE: Algorithm failed with code {}".format(self.status))
            return
        else:
            self.transects_valid_count = 0
            sum_valid = 0
            count_all = 0
            sum_all = 0
            windspeed = sqrt(self.u * self.u + self.v * self.v) # (m/s)

            direction_correction = cos(pi * (self.downwind_box_azimuth_adjusted - self.downwind_box_azimuth_initial) / 180.0)
            if (not(logger is None)): logger.info("Direction correction {}".format(direction_correction))
            for t in self.transects:
                sum_all += t.enhkg * direction_correction
                count_all += 1
                if (t.isvalid()):
                    sum_valid += t.enhkg * direction_correction
                    self.transects_valid_count += 1

            if (self.transects_valid_count < 3):
                # Calculate q in Gg/hour 
                average_all = sum_all/count_all # (kg/m)
                q_all= windspeed * average_all * 1E-6 * 3600.0
                if (not(logger is None)): logger.info("All transect count: {} Average enhancement over transect {} (kg/m)".format(count_all, average_all))
                if (not(logger is None)): logger.info("Windspeed: {} Q {} (Gg/hour)".format(windspeed, q_all))
                if (not(logger is None)): logger.info("REJECTED: Count of valid transects: {} is less than required minimum of 3.".format(self.transects_valid_count))

                if self.downwind_box_count == 0:
                    self.status |= Status.FAILURE_NODATA_DOWNWINDBOX.value
                    if (not(logger is None)): logger.info("REJECTED: No data in downwindbox. Count: {}".format(self.transects_valid_count))

                if windspeed < 2: 
                    self.status |= Status.FAILURE_WIND_BELOW2.value
                    if (not(logger is None)): logger.info("REJECTED: Windspeed below 2 m/s. Speed: {}".format(windspeed))


                self.status |= Status.FAILURE_TRANSECTS_BELOW_3.value
                if (not(logger is None)): logger.info("REJECTED: Insufficient number of valid transects. Count: {}".format(self.transects_valid_count))
            
            else:
                if windspeed < 2: 
                    self.status |= Status.FAILURE_WIND_BELOW2.value
                    if (not(logger is None)): logger.info("REJECTED: Windspeed below 2 m/s. Speed: {}".format(windspeed))

                if self.downwind_box_count == 0:
                    self.status |= Status.FAILURE_NODATA_DOWNWINDBOX.value
                    if (not(logger is None)): logger.info("REJECTED: No data in downwindbox. Count: {}".format(self.transects_valid_count))

                if self.transects_valid_count < 2: 
                    self.status |= Status.FAILURE_TRANSECTS_BELOW_3.value
                    if (not(logger is None)): logger.info("REJECTED: Insufficient number of valid transects. Count: {}".format(self.transects_valid_count))

                average_valid = sum_valid/self.transects_valid_count # (kg/m)
                # Calculate q in Gg/hour 
                self.q_valid_hour = windspeed * average_valid * 1E-3 * 3600.0 # t/hour
                self.q_valid_year = self.q_valid_hour * 365 * 24 * 1E-3 # Gg/year
                if (not(logger is None)): logger.info("Valid transect count: {}".format(self.transects_valid_count))
                if (not(logger is None)): logger.info("Average enhancement over transect: {:.3f} (kg/m)".format(average_valid))
                if (not(logger is None)): logger.info("Windspeed: {} (m/s)".format(windspeed))
                if (not(logger is None)): logger.info("Q: {:.3f} (t/hour) = {:.3f} (Gg/year)".format(self.q_valid_hour, self.q_valid_year))

                if self.status == Status.PROCESSING.value:
                    self.status = Status.SUCCESS.value
                    if (not(logger is None)): logger.info("SUCCESS")
                else:
                    if self.q_valid_hour is None:
                        print("Orbit: {} FAILURE: {} Q: NA (t/hour) = NA (Gg/year)".format(self.tropomi.orbit, self.status))
                    else:
                        print("Orbit: {} FAILURE: {} Q: {:.3f} (t/hour) = {:.3f} (Gg/year)".format(self.tropomi.orbit, self.status, self.q_valid_hour, self.q_valid_year))
  
    def calculate_transects(self):
        '''
        From Sadavarte 2021: 
        We define 15 equally spaced transects between the source and the end of the rectangular mask for calculating the source rates.
        For each transect populate the following attributes
        enhkg - Enhancement (kg) - Sum over pixels not masked (pixel value - background) x (length of intersect)
        len_mask - Length of intersection with pixels masked (m)
        len_values - Length of intersection with pixels not masked (m)
        pix_mask - Count of intersection with pixels masked 
        pix_values - Count of intersection with pixels not masked
        '''

        if not (self.status == Status.PROCESSING.value): return
        if (self.downwind_box_half_width == None): self._calculate_downwindbox_widthadjustment()

        self.transects: list[Transect] = []
        d14 = self.downwind_box_length / 14.0
        # Repeat length adjustment procedure starting from 4th line
        for i in range(3,15):
            [[pm_lon, pm_lat, _]] = Algorithm_CSF.g.direct(points=[self.source.xy.x, self.source.xy.y], azimuths=self.downwind_box_azimuth_adjusted, distances=d14 * i)
            pm = [pm_lon, pm_lat]
            [[ps_lon, ps_lat, _]] = Algorithm_CSF.g.direct(points=pm, azimuths=self.downwind_box_azimuth_adjusted - 90, distances=self.downwind_box_half_width)
            [[pe_lon, pe_lat, _]] = Algorithm_CSF.g.direct(points=pm, azimuths=self.downwind_box_azimuth_adjusted + 90, distances=self.downwind_box_half_width)
            self.transects.append(Transect(ps_lon, ps_lat,pe_lon, pe_lat))
            
        [minscan, minpixel, maxscan, maxpixel, lonbox, latbox, _, _, scan_downwindbox]  = self.tropomi.narrow_to_domain(Geometry.INTERSECTS, self.downwind_box)
        scans = maxscan - minscan
        pixels = maxpixel - minpixel
        pressure = self.tropomi.get_variable_data("surface_pressure")[0, minscan:maxscan, minpixel:maxpixel]
        dry_air = self.tropomi.get_variable_data("dry_air_subcolumns")[0,minscan:maxscan, minpixel:maxpixel,:]
        h2o = self.tropomi.get_variable_data("water_total_column")[0, minscan:maxscan, minpixel:maxpixel]
        self.downwind_box_minimumvalue = ma.min(scan_downwindbox)
        if self.config.Algorithm_CSF_transect_background == "minimum": self.background = self.downwind_box_minimumvalue

        # For each pixel in this domain
        for s in range(scans):
            for p in range(pixels):
                pixel = Polygon([
                    (lonbox[s,p,0], latbox[s,p,0]),
                    (lonbox[s,p,1], latbox[s,p,1]),
                    (lonbox[s,p,2], latbox[s,p,2]),
                    (lonbox[s,p,3], latbox[s,p,3])
                ])
                if not(self.downwind_box.intersects(pixel)): continue
                isvalid = False
                isvalidpositive = False
                for t in self.transects:
                    line = t.line
                    if pixel.intersects(line):
                        int_l = pixel.intersection(line).coords
                        ps = [int_l[0]]
                        pe = [int_l[1]]
                        [[dist, _, _]] = Algorithm_CSF.g.inverse(ps, pe)
                        
                        if scan_downwindbox.mask[s,p]:
                            t.len_mask += dist
                            t.pix_mask += 1
                        else:
                            isvalid = True
                            v = scan_downwindbox[s,p] - self.background
                            t.len_values += dist
                            t.pix_values += 1
                            # Add only pixels with positive enhancement
                            if (0 < v): 
                                t.enhkg += dist *  convert_column_ppb_with_water(v, average(dry_air[s,p,:]), h2o[s,p], pressure[s,p])
                                t.len_values_positive += dist
                                t.pix_values_positive += 1
                                isvalidpositive = True
                if isvalid: self.transects_valid_pixels_count += 1
                if isvalidpositive: self.transects_valid_positive_pixels_count += 1

        for t in self.transects:
            if (not(logger is None)): logger.info("{}".format(t))

        return

    def count_source_in_downindbox(self) -> int:
        ''' Calculate number of sources in the downwind box if exists.'''
        cnt = 0
        if self.downwind_box is None: return cnt
        for s in Source.Sources:
            src = create_source(s)
            if src == self.source: continue
            if self.downwind_box.contains(src.xy): cnt += 1
        return cnt
    
    def count_source_in_upwindbox(self) -> int:
        ''' Calculate number of sources in the upwind box if used.'''
        cnt = 0
        if self.upwind_box is None: return cnt
        for s in Source.Sources:
            src = create_source(s)            
            if src == self.source: continue
            if self.upwind_box.contains(src.xy): cnt += 1
        return cnt

    def list_calculation_Properties(self):
        '''
        List calculation properties for data analysis.
        Entry point into this algorithm
        '''
        print("Orbit: {} Processor: {} Date: {}".format(self.tropomi.orbit, self.tropomi.processor_version, datetime.strftime(self.tropomi_source_date,"%Y-%m-%d %H:%M")))
        print("Wind: u: {} (m/s) v: {} (m/s) speed: {} (m/s)".format(self.u, self.v, sqrt(self.u * self.u + self.v * self.v)))

        print("Configuration")
        print("\tBackground geometry: {}".format(self.config.Algorithm_CSF_background_geometry))
        print("\tBackground minimum pixel count: {}".format(self.config.Algorithm_CSF_background_minimum_count))
        print("\tDownwind box geometry: {}".format(self.config.Algorithm_CSF_downwindbox_geometry))
        print("\tDownwind box mask: {}".format(self.config.Algorithm_CSF_downwindbox_mask))
        print("\tDownwind box enhancement delta: {}".format(self.config.Algorithm_CSF_downwindbox_enhancement_delta))
        print("\tTransect geometry: {}".format(self.config.Algorithm_CSF_transect_background))

        print("Background")
        print("\tGeometry: {}".format(self.config.Algorithm_CSF_background_geometry.name))
        print("\tUpwind box: Count: {} Upwind Average: {}".format(self.upwind_box_background_count, self.upwind_box_background_average))
        print("\tDomain median: {} ".format(self.domain_median))
        if not self.config.Algorithm_CSF_transect_background_value is None: print("\tManual background: {}".format(self.config.Algorithm_CSF_transect_background_value))
        else: print("\tManual background: Not used")
        print("\tBackground enhancement: {}".format(self.background_enhancement))
        print("\tBackground downwind box: {}".format(self.background))

        print("Downwind box")
        print("\tGeometry: {} Mask: {}".format(self.config.Algorithm_CSF_downwindbox_geometry.name, self.config.Algorithm_CSF_downwindbox_mask.name))
        print("\tRotation: {}".format(self.downwind_box_azimuth_adjusted - self.downwind_box_azimuth_initial))
        print("\tLength: {} (m) Width: {} (m)".format(self.downwind_box_length, self.downwind_box_half_width * 2))
        
        print("Transects")
        print("\tValid pixels count: {} Valid positive pixel count: {}".format(self.transects_valid_pixels_count, self.transects_valid_positive_pixels_count))
        print("\tValid transect count: {} Minimum: {} (ppb)".format(self.transects_valid_count, self.downwind_box_minimumvalue))
        for t in self.transects:
            print("\t{}".format(t))

        if not (self.q_valid_hour is None):
            print("Emission rate: {} (t/hour) {} (Gg/year)".format(self.q_valid_hour, self.q_valid_year))
        else:
            print("Emission rate has not been calculated")
        for s in Status:
            if s.value & self.status or s.value == self.status: print("Status: {}".format(s.name))

    @classmethod
    def run(cls, config: Config, source: Source, orbit:str, processor: str, wind: tuple[float]):
        '''
        Execute algorithm from begining or from the last completed step. Entry point into algorithm

        Parameters
        ----------
        config: Config
            Run configuration
        source: Source
            Source of emissions
        orbit: str
            orbit of emissions
        processort: str
            processor for orbit
        wind: tuple
            u,v components of boundary layer pressue averaged wind

        Returns
        -------
        cls: Sadavarte2021
            Initialized object
        '''
        from os import path
        filename_end = Algorithm_CSF.get_picklename(config, source.case_name, orbit, processor)

        # Restart calculation at the last succesfull place
        if (path.exists(filename_end)):
            with open(filename_end, mode="rb") as f:
                cls: Algorithm_CSF = load(f)
                f.close()

            # It is the same configuration do nothing
            if cls.config == config:
                if not (logger is None): logger.info("Previous execution was found. ")
                return cls
            
            # Configuration changed so need to clean up pickle files
            if path.exists(filename_end): remove(filename_end)               

        tropomi = TROPOMI_for_orbit(orbit, processor)
        cls = Algorithm_CSF(config, source, wind, tropomi)
        cls._calculate_background()
        cls._calculate_downwindbox_directionadjustment()
        cls._calculate_downwindbox_lengthadjustment()
        cls._calculate_downwindbox_widthadjustment()
        cls._calculate_transects()
        cls.calculate_sourcerate()

        if cls.status == Status.SUCCESS.value: 
            print("Orbit: {} SUCCESS: Q: {:.3f} (t/hour) = {:.3f} (Gg/year)".format(cls.tropomi.orbit, cls.q_valid_hour, cls.q_valid_year))
        else:
            if cls.q_valid_hour is None:
                if cls.status & Status.FAILURE_WIND_BELOW2.value: print("Orbit: {} FAILURE: Wind below 2(m/s).".format(cls.tropomi.orbit))
                if cls.status & Status.FAILURE_BACKGROUND_VALUES_HIGH.value: print("Orbit: {} FAILURE: Background box contains no data.".format(cls.tropomi.orbit))
                if cls.status & Status.FAILURE_NODATA_DOWNWINDBOX.value: print("Orbit: {} FAILURE: Downwind box contains no data.".format(cls.tropomi.orbit))
                if cls.status & Status.FAILURE_TRANSECTS_BELOW_3.value: print("Orbit: {} FAILURE: Fewer then 3 valid transects.".format(cls.tropomi.orbit))
                if cls.status & Status.FAILURE_NODATA_INDOMAIN.value: print("Orbit: {} FAILURE: No data in domain.".format(cls.tropomi.orbit))
                if cls.status & Status.FAILURE_TRANSECTS_PIXEL_COUNT.value: print("Orbit: {} FAILURE: Number of valid pixels intersecting transects below threshold.".format(cls.tropomi.orbit))
                if cls.status & Status.FAILURE_ROTATION_0_NO_PIXELS.value: print("Orbit: {} FAILURE: Rotation 0 does not contain any pixels.".format(cls.tropomi.orbit))
            else:
                if cls.status & Status.FAILURE_WIND_BELOW2.value: print("Orbit: {} FAILURE: Wind below 2(m/s). Q: {:.3f} (t/hour) = {:.3f} (Gg/year)".format(cls.tropomi.orbit, cls.q_valid_hour, cls.q_valid_year))
                else: print("Orbit: {} FAILURE: {} Q: {:.3f} (t/hour) = {:.3f} (Gg/year)".format(cls.tropomi.orbit, cls.status, cls.q_valid_hour, cls.q_valid_year))
            
        with open(cls.picklename, mode="wb") as f:
            dump(cls, f)
            f.close()
        return cls
 
    @classmethod
    def get_picklename(cls, config: Config, source: str, orbit: str, processor: str) -> str:
        '''
        Construct path for a pickle file
        '''
        temp = "Algorithm_CSF_{}_{}_{}{}{}.pkl".format(orbit, processor, config.Algorithm_CSF_background_geometry.value, config.Algorithm_CSF_downwindbox_geometry.value, config.Algorithm_CSF_downwindbox_mask.value)
        return path.join(Config.CSF_folder, source, temp )

def _handler_chart_all(args: Namespace):
    filename_end = path.join(Config.CSF_folder, args.source, args.orbit, args.processor, "Algorithm_CSF_{}_{}.pkl".format(args.source, args.orbit))
    if not path.exists(filename_end):
        print("Error. No model run {}. Please run algorithm first".format(filename_end))
        exit()

    with open(filename_end, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sad._chart_combined(None)
    sad.chart_lonlat_valid(None)
    sad.chart_lonlat_positive(None)
    sad._chart_downwindbox_imagecoords_valid(None)
    sad._chart_downwindbox_imagecoords_positive(None)
    sad.chart_transects(None)

def _handler_chart_combined(args: Namespace):
    config = Config()
    filename_end = Algorithm_CSF.get_picklename(config, args.source, args.orbit, args.processor)
    if not path.exists(filename_end):
        print("Error. No model run {}. Please run algorithm first".format(filename_end))
        exit()

    with open(filename_end, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sad._chart_combined(None)

def _handler_chart_positive(args: Namespace):
    config = Config()
    filename_end =  Algorithm_CSF.get_picklename(config, args.source, args.orbit, args.processor)
    if not path.exists(filename_end):
        print("Error. No model run {}. Please run algorithm first".format(filename_end))
        exit()

    with open(filename_end, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sad.chart_lonlat_positive(None)
    return

def _handler_chart_positive_image(args: Namespace):
    config = Config()
    filename_end =  Algorithm_CSF.get_picklename(config, args.source, args.orbit, args.processor)
    if not path.exists(filename_end):
        print("Error. No model run {}. Please run algorithm first".format(filename_end))
        exit()

    with open(filename_end, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sad._chart_downwindbox_imagecoords_positive(None)
    return

def _handler_chart_transects(args: Namespace):
    config = Config()
    filename_end =  Algorithm_CSF.get_picklename(config, args.source, args.orbit, args.processor)
    if not path.exists(filename_end):
        print("Error. No model run {}. Please run algorithm first".format(filename_end))
        exit()

    with open(filename_end, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sad.chart_transects(None)

def _handler_chart_valid(args: Namespace):
    config = Config()
    filename_end =  Algorithm_CSF.get_picklename(config, args.source, args.orbit, args.processor)
    if not path.exists(filename_end):
        print("Error. No model run {}. Please run algorithm first".format(filename_end))
        exit()

    with open(filename_end, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sad.chart_lonlat_valid(None)
    return

def _handler_chart_valid_image(args: Namespace):
    config = Config()
    filename_end =  Algorithm_CSF.get_picklename(config, args.source, args.orbit, args.processor)
    if not path.exists(filename_end):
        print("Error. No model run {}. Please run algorithm first".format(filename_end))
        exit()

    with open(filename_end, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sad._chart_downwindbox_imagecoords_valid(None)
    return

def _handler_clean(args: Namespace):
    '''
    Handler for clean menu
    '''
    partname = "Algorithm_CSF_{}_{}".format(args.source, args.orbit)
    config = Config()
    filename_end =  Algorithm_CSF.get_picklename(config, args.source, args.orbit, args.processor)
    filename_log = path.join("Log", partname + ".log")

    if path.exists(filename_log): remove(filename_log)
    if path.exists(filename_end): remove(filename_end)

def _handler_clean_all(args: Namespace):
    '''
    Handler for clean menu
    '''
    dir_list = listdir(path.join(Config.CSF_folder, args.source))
    for orbit in dir_list:        
        filename_log = path.join("Log", "Algorithm_CSF_{}_{}.log".format(args.source, orbit))
        if path.exists(filename_log): 
            remove(filename_log)
            if not logger is None: logger.info("Removed: {}".format(filename_log))

        cache = path.join(Config.CSF_folder, args.source)
        if not path.exists(cache): continue
        if not path.isdir(cache): continue

        for p in Config.TROPOMI_Processors:
            cache_list = listdir(cache)
            processor = "_"+p+"_"
            for pkl in cache_list:
                if pkl.startswith("Algorithm_CSF_") and processor in pkl and pkl.endswith(".pkl"):
                    filename_end = path.join(cache, pkl)        
                    if path.exists(filename_end): 
                        remove(filename_end)
                        if not logger is None: logger.info("Removed: {}".format(filename_end))

def _handler_list(args: Namespace):
    '''
    Handler for list run menu

    Parameters
    ----------
    args.source:str
        Name of the source (mandatory)
    args.orbit:str
        Orbit (mandatory)
    '''
    filename_end = Algorithm_CSF.get_picklename(Config(), args.source, args.orbit, args.processor)

    if (path.exists(filename_end)):
        with open(filename_end, mode="rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()
        sad.list_calculation_Properties()
    else:
        print("ERROR: Unable to find file {}".format(filename_end))

    return None   

def _handler_run(args: Namespace):
    '''
    Handler for run menu
    
    Parameters
    ----------
    args.source (str) :
        Name of the source (mandatory)
    args.orbit:str
        Orbit (mandatory)
    args.processor:str
        Processor (mandatory)
    args.transect_background_method:str
        Method of transect background calculations allowing manual override of algorithm
    args.transect_background_manual_value:str
        Background value mandatory when bacgroundtype is set to manual
    '''
    from ERA5 import ERA5, ERA5_PressureLevels, ERA5_SingleLevel
    
    config = Config()
    if not args.transect_background_method is None: 
        config.Algorithm_CSF_transect_background = args.transect_background_method
        if args.transect_background_method == "manual":
            if args.transect_background_manual_value is None:
                message = "Error: When background type is set to 'manual'. Background value is mandatory"
                if not(logger is None): logger.error(message)
                print(message)
                return
            else:
                config.Algorithm_CSF_transect_background_value = float(args.transect_background_manual_value)
    
    if not args.background_geometry is None: 
        if args.background_geometry == "contains": config.Algorithm_CSF_background_geometry = Geometry.CONTAINS
        if args.background_geometry == "center": config.Algorithm_CSF_background_geometry = Geometry.CENTER
        if args.background_geometry == "intersects": config.Algorithm_CSF_background_geometry = Geometry.INTERSECTS

    if not args.downwindbox_geometry is None:
        if args.downwindbox_geometry == "contains": config.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
        if args.downwindbox_geometry == "center": config.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
        if args.downwindbox_geometry == "intersects": config.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS

    if not args.downwindbox_mask is None:
        if args.downwindbox_mask == "none": config.Algorithm_CSF_downwindbox_mask = Mask.NONE
        if args.downwindbox_mask == "negative": config.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE

    source = create_source(args.source)
    Config.create_log("Algorithm_CSF_{}_{}".format(args.source, args.orbit) + ".log")

    tropomi = TROPOMI_for_orbit(args.orbit, args.processor)

    Config.create_directory_structure_CSF(args.source)
    [_,_,dt, _] = tropomi.get_pixel_for_source(source.xy)
    grid_time = ERA5.calculate_ERA_GridTime(dt)
    
    path_ERA5_SL = ERA5_SingleLevel.get_filepath("CSF", args.source, args.orbit)
    if (not(path.exists(path_ERA5_SL))):
        message = "Required file: {} does not exist. Please download this file first using ERA5 program.".format(path_ERA5_SL)
        if not logger is None: logger.error(message)
        print(message)
        return

    era5_sl = ERA5_SingleLevel.from_netCDF(path_ERA5_SL)
    
    path_ERA5_PL = ERA5_PressureLevels.get_filepath("CSF", source.case_name, args.orbit)
    if not(path.exists(path_ERA5_PL)) :
        message = "Required file: {} does not exist. Please download this file first using ERA5 program.".format(path_ERA5_PL)
        if not logger is None: logger.error(message)
        print(message)
        return
    
    era5_p = ERA5_PressureLevels.from_netCDF(path_ERA5_PL)

    era5_p.convert_timestamp_to_index(grid_time)
    era5_p.convert_lat_to_index(source.gridpoint.y)
    era5_p.convert_lon_to_index(source.gridpoint.x)

    # Pressure, geopotential, u and v as 1-d arrays over pressure levels
    p = era5_p.pressure
    z = era5_p.variables["z"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]
    u = era5_p.variables["u"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]
    v = era5_p.variables["v"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]


    era5_sl.convert_timestamp_to_index(grid_time)
    era5_sl.convert_lat_to_index(source.gridpoint.y)
    era5_sl.convert_lon_to_index(source.gridpoint.x)
    blh = era5_sl.get("blh", era5_sl.time_index, era5_sl.lat_index, era5_sl.lon_index) 
    sp = era5_sl.get("sp", era5_sl.time_index, era5_sl.lat_index, era5_sl.lon_index) / 100.0
    
    wind = Meteorology.calculate_ERA5_PressureAveragedWind(p, z, u, v, blh, sp)
    _ = Algorithm_CSF.run(config, source, args.orbit, args.processor, wind)

if __name__ == '__main__':
    '''
    Entry point into this module.
    Parses command line arguments and call handlers for different menu options
    To find command options execute: python Algorithm_CSF.py -h
    '''
    help_source = "Source of emissions"
    help_orbit = "Orbit number as a five digit string e.g. 04579"
    help_mask = "Masking option."
    
    parser = ArgumentParser(prog="Algorithm_CSF", description="Implementation of Cross Sectional FLux method following Saadavarte et all. 2021")
    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_chart = subparsers.add_parser("chart", help="Charts for algorithm.")
    subparsers_chart = parser_chart.add_subparsers(help="subcommand help", required=True)

    parser_chart_all = subparsers_chart.add_parser("all", help="All")
    parser_chart_all.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_chart_all.add_argument("orbit", type=str, help=help_orbit)
    parser_chart_all.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart_all.set_defaults(func=_handler_chart_all)

    parser_chart_combined = subparsers_chart.add_parser("combined", help="Positive and Valid together in geographical coordinates")
    parser_chart_combined.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_chart_combined.add_argument("orbit", type=str, help=help_orbit)
    parser_chart_combined.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart_combined.set_defaults(func=_handler_chart_combined)

    parser_chart_positive = subparsers_chart.add_parser("positive", help="Positive only in geographical coordinates")
    parser_chart_positive.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_chart_positive.add_argument("orbit", type=str, help=help_orbit)
    parser_chart_positive.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart_positive.set_defaults(func=_handler_chart_positive)

    parser_chart_positive_image = subparsers_chart.add_parser("positive_image", help="Positive only in image coordinates")
    parser_chart_positive_image.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_chart_positive_image.add_argument("orbit", type=str, help=help_orbit)
    parser_chart_positive_image.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart_positive_image.set_defaults(func=_handler_chart_positive_image)

    parser_chart_transects = subparsers_chart.add_parser("transects", help="Emission rate for valid transect")
    parser_chart_transects.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_chart_transects.add_argument("orbit", type=str, help=help_orbit)
    parser_chart_transects.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart_transects.set_defaults(func=_handler_chart_transects)

    parser_chart_valid = subparsers_chart.add_parser("valid", help="Valid in geographical coordinates")
    parser_chart_valid.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_chart_valid.add_argument("orbit", type=str, help=help_orbit)
    parser_chart_valid.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart_valid.set_defaults(func=_handler_chart_valid)

    parser_chart_valid_image = subparsers_chart.add_parser("valid_image", help="Valid in image coordinates")
    parser_chart_valid_image.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_chart_valid_image.add_argument("orbit", type=str, help=help_orbit)
    parser_chart_valid_image.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_chart_valid_image.set_defaults(func=_handler_chart_valid_image)

    parser_clean = subparsers.add_parser("clean", help="Clean cache and log for specified source and orbit.")
    parser_clean.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_clean.add_argument("orbit", type=str, help=help_orbit)
    parser_clean.set_defaults(func=_handler_clean)

    parser_clean_all = subparsers.add_parser("cleanall", help="Clean cache and log for specified source and all orbits.")
    parser_clean_all.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_clean_all.set_defaults(func=_handler_clean_all)

    parser_list = subparsers.add_parser("list", help="List algorithm details for a run.")    
    parser_list.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_list.add_argument("orbit", type=str, help=help_orbit)
    parser_list.add_argument("processor", type=str, help="Processor version as six digits")
    parser_list.set_defaults(func=_handler_list)

    parser_run = subparsers.add_parser("run", help="run algorithm for source and orbit with geometry and masking")
    parser_run.add_argument("source", type=str, choices=Source.Sources, help=help_source)
    parser_run.add_argument("orbit", type=str, help=help_orbit)
    parser_run.add_argument("processor", type=str, help="Processor version as six digits")
    parser_run.add_argument("-bg", "--background_geometry", type=str, choices=["contains", "center", "intersects"], help="Overwrite default geometry used for background calculation.")
    parser_run.add_argument("-dg", "--downwindbox_geometry", type=str, choices=["contains", "center", "intersects"], help="Overwrite default geometry used for downwindbox calculation.")
    parser_run.add_argument("-dm", "--downwindbox_mask", type=str, choices=["none", "negative"], help="Overwrite default masking used for downwindbox calculation.")
    parser_run.add_argument("-tbm", "--transect_background_method", type=str, choices=["algorithm", "domain", "manual", "upwind", "minimum"], help="Overwrite default choice of background used for transect calculation")
    parser_run.add_argument("-tbv", "--transect_background_manual_value", type=float, default=None, help="Value used as background for backgroundtype=manual" )
    parser_run.set_defaults(func=_handler_run)

    args = parser.parse_args()    
    argumentsvalid = True

    if ("orbit" in args):
        if not(len(args.orbit) == 5) or not(args.orbit.isdigit()):
            print ("Orbit must be a five digit number: {}".format(args.orbit))
            argumentsvalid = False

    if ("processor" in args):
        if not(len(args.processor) == 6) or not(args.processor.isdigit()):
            print ("Processor must be a 6 digit number: {}".format(args.processor))
            argumentsvalid = False
    
    if argumentsvalid: args.func(args)
