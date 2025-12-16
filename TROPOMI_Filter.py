from argparse import ArgumentParser, Namespace
from Config import Config
from Copernicus_Catalogue import Copernicus_Catalogue
from datetime import datetime, timedelta, timezone
from Geometry import Geometry
import logging
from math import fabs, sqrt
from numpy import ma
from os import listdir, path, remove
from pickle import dump, load
from shapely import Point, Polygon
from TROPOMI import TROPOMI

logger = logging.getLogger(__name__)

geometry = Geometry.CONTAINS

class TROPOMI_Filter:
    '''
    Progress and results of filter operation as specified in Sadavarte 2021
    '''
    pickle = path.join(Config.Paper_folder, "TROPOMI_Filter_HailCreek.pkl")
    
    def __init__(self, step:str, config: Config, date_first_inclusive:datetime, date_last_exclusive:datetime):
        '''
        Intialize filter. The first filter to run is contains

        Parameters
        ----------
        step:str
            Filter step. First step is "contains", second "qa", last one "r". Corresponding to algorithm from Sadavarte
        domain:Polygon
            domain for filter
        orbit_first_inclusive:int
            First orbit to be filtered
        orbit_last_exclusive:int
            The first orbit not to be filtered
        '''
        self.config = config
        self.step = step

        # These fields are used by "qa" and "r" filters
        # orbit_first_inclusive should be set by contains filter
        self.orbit_first_inclusive:int = None
        # orbit_last_exclusive should be set by contains filter
        self.orbit_last_exclusive:int = None

        self.orbitprocessed_qa: int = None
        self.success_qa: list[int] = []

        self.orbitprocessed_r: int = None
        self.success_r: list[int] = []

        # These fields are used by "contains" filter that works on description files using contains filter
        self.date_first_inclusive:datetime = date_first_inclusive
        self.date_last_exclusive:datetime = date_last_exclusive

        self.dateprocessed: datetime = None
        self.success_contains: list[int] = []
        self.success_contains_s3_key: list[str] = []
        # Current filter results
        self.result: list[int] = []
    
    def completed(self, step: str):
        '''
        True if run through all days or all orbits
        '''
        if step == "contains":
            if self.dateprocessed is None: return False
            td1 = timedelta(days=1)
            return self.date_last_exclusive <= self.dateprocessed + td1
        if step == "qa":
            if self.orbitprocessed_qa is None: return False
            return self.orbit_last_exclusive <= self.orbitprocessed_qa + 1
        if step == "r":
            if self.orbitprocessed_r is None: return False
            return self.orbit_last_exclusive <= self.orbitprocessed_r + 1
    
    def _filter_1_contains(self, username: str, password:str):
        '''
        Select TROPOMI files containg domain and processed by speficin versions of processor
        '''
        catalogue = Copernicus_Catalogue(username, password, self.config)
        date = self.date_first_inclusive if self.dateprocessed is None else self.dateprocessed

        while date < self.date_last_exclusive:
            # This returns a list holding AWS s3 bucket keys
            # /Sentinel-5P/TROPOMI/L2__CH4___/2018/12/08/S5P_RPRO_L2__CH4____20181208T022821_20181208T040950_05969_03_020400_20221117T000403.nc
            # /Sentinel-5P/TROPOMI/L2__CH4___/2018/12/08/S5P_OFFL_L2__CH4____20181208T022821_20181208T040950_05969_01_010202_20181214T041145.nc
            s3files = catalogue.get_s3_keys(date.strftime("%Y%m%d"))
            for d in s3files:
                f = d[-86:]
                iorbit=TROPOMI.get_orbit_int(f)
                processor = TROPOMI.get_processor_version(f)
                
                # Skip files that we are not interested in
                if not processor in Config.TROPOMI_Processors: continue

                if not d in self.success_contains_s3_key:
                    if not logger is None: logger.info("Adding orbit:{} s3 key: {}".format(iorbit, d))
                    if not iorbit in self.success_contains: self.success_contains.append(iorbit)
                    self.success_contains_s3_key.append(d)
                
            if self.orbit_first_inclusive is None or iorbit < self.orbit_first_inclusive : self.orbit_first_inclusive = iorbit
            if self.orbit_last_exclusive is None or self.orbit_last_exclusive <= iorbit : self.orbit_last_exclusive = iorbit + 1

            self.dateprocessed = date
            self.result = self.success_contains.copy()
            with open(self.pickle, mode="wb") as f:
                dump(self, f)
                f.close()
            
            date = date + timedelta(days = 1) 

    def _filter_2_quality(self):
        '''
        Find all scans from April 30 2018 to December 31 2019 containing at least 500 pixels in domain 20-24S 146-150 E
        From Sadavarte 2021: For emission quantification from TROPOMI-detected plumes, orbits from 2018 and 2019 were 
        screened with >500 individual observation pixels in the domain of 20° - 24°S and 146° - 150°E
        '''

        if self.step == "contains" and not self.completed("contains"):
            print("The last orbit processed by 'Contains' filter was {} and has not completed. Filter needs to run until {}.".format(self.dateprocessed, self.date_last_exclusive))
            exit()

        # This filter contains names of description files their names are identical to TROPOMI but end with .xml
        self.step = "qa"
        orbits = self.success_contains.copy()

        processors = Config.TROPOMI_Processors.copy()
        processors.sort(reverse=True)
        if not(logger is None): logger.info("Start: {} End: {}".format(ma.min(orbits), ma.max(orbits)))

        ncfiles = self.find_highest_configured_procesor(processors)
        
        for i in range(self.orbit_first_inclusive, self.orbit_last_exclusive):
            # Go to the last processed orbit
            if not (self.orbitprocessed_qa is None) and i <= self.orbitprocessed_qa: continue

            self.orbitprocessed_qa = i
            if not i in orbits: 
                logger.info("OMITTED: NC orbit: {} is not included in the initial set".format(i))
                with open(self.pickle, mode="wb") as f:
                    dump(self, f)
                    f.close()
                continue

            if not (i in ncfiles): 
                logger.info("OMITTED: NC file for orbit: {} is not available".format(i))
                with open(self.pickle, mode="wb") as f:
                    dump(self, f)
                    f.close()
                continue

            f = ncfiles[i]
            processor = TROPOMI.get_processor_version(f)
            filename = path.join(Config.TROPOMI_folder, processor, f)
            tropomi = TROPOMI(filename)
            [minscan, minpixel, maxscan, maxpixel, ln_box_arr, lt_box_arr, lons, lats, scan] = tropomi.narrow_to_domain(geometry, self.config.Algorithm_CSF_domain)
            if minscan is None or minpixel is None or maxscan is None or maxpixel is None:
                logger.info("No pixels found for orbit: {}".format(i))
            else:
                scans = maxscan - minscan
                pixels = maxpixel - minpixel
                if (ma.count(scan) < 500):
                    logger.info("Fewer then 500 pixels in domain found for orbit: {}".format(i))
                else:
                    qa_value = tropomi.get_variable_data("qa_value")[0, minscan:maxscan , minpixel:maxpixel]
                    status = "REJECTED"

                    '''
                        From Sadavarte 2021: The TROPOMI XCH4 measurements used in this analysis were screened for cloud-free coverage and 
                        low aerosol content using the quality flag provided in the data products (we use qa = 1).
                        Data quality qa = 1 signifies XCH4 is filtered for solar zenith angle (<70°), viewing zenith angle (<60°), 
                        smooth topography (1 standard deviation surface elevation variability <80 m within a 5 km radius), and 
                        low aerosol load (aerosol optical thickness <0.3 in the NIR band). 
                        The TROPOMI data was corrected for XCH4 variations due to surface elevation by adding 7 ppb per km surface elevation with
                        respect to the mean sea level.25 TROPOMI XCH4 data show artificial stripes in the flight direction, most probably due to
                        swath position-dependent calibration inaccuracies, which were corrected by applying a fixed mask destriping approach to the L2
                        data developed for the TROPOMI XCO retrieval. 
                        For emission quantification from TROPOMI-detected plumes, orbits from 2018 and 2019 were screened with >500 individual 
                        observation pixels in the domain of 20°-24°S and 146°-150°E (Figure 1a). 
                        To ensure that emission quantifications are not influenced by systematic surface albedo or aerosol bias, we reject orbits that 
                        show a high correlation (|R| > 0.5) of XCH4 with surface albedo or aerosol optical thickness.
                    '''
                    counter = 0
                    for s in range(scans):
                        for p in range(pixels):

                            if qa_value[s,p] < 1.0: continue

                            if (geometry == Geometry.INTERSECTS):
                                # See file:///Z:/Sentinel5/Documentation/Sentinel-5P-Level-2-Input-Output-Data-Definition.pdf page 141
                                ul = (ln_box_arr[s,p,3], lt_box_arr[s,p,3])
                                bl = (ln_box_arr[s,p,0], lt_box_arr[s,p,0])
                                br = (ln_box_arr[s,p,1], lt_box_arr[s,p,1])
                                ur = (ln_box_arr[s,p,2], lt_box_arr[s,p,2])
                                pixelbox = Polygon([ul,bl,br,ur])
                                if not self.config.Algorithm_CSF_domain.intersects(pixelbox): continue
                            elif (geometry == Geometry.CONTAINS) :
                                ul = (ln_box_arr[s,p,3], lt_box_arr[s,p,3])
                                bl = (ln_box_arr[s,p,0], lt_box_arr[s,p,0])
                                br = (ln_box_arr[s,p,1], lt_box_arr[s,p,1])
                                ur = (ln_box_arr[s,p,2], lt_box_arr[s,p,2])
                                pixelbox = Polygon([ul,bl,br,ur])
                                if not self.config.Algorithm_CSF_domain.contains(pixelbox): continue
                            elif (geometry == Geometry.CENTER):
                                center = Point(lons[s,p], lats[s,p])
                                if not self.config.Algorithm_CSF_domain.contains(center): continue
                            
                            counter += 1
                    
                    if 500 < counter: 
                        status = "ACCEPTED"
                        self.success_qa.append(i)
                    if not(logger is None):
                        message = "Orbit={}, CONTAINS: {}, Count={}".format(i, status, counter)
                        logger.info(message)

            self.result = self.success_qa.copy()
            with open(self.pickle, mode="wb") as f:
                dump(self, f)
                f.close()
        pass
 
    def _filter_3_correlation(self):
        '''
        To ensure that emission quantifications are not influenced by systematic surface albedo or aerosol
        bias, we reject orbits that show a high correlation (|R| > 0.5) of
        XCH4 with surface albedo or aerosol optical thickness.
        '''
        
        if self.step == "qa" and not self.completed("qa"):
            print("The last orbit processed by 'qa' filter was {}. Filter needs to run until {}.".format(self.orbitprocessed_qa, self.orbit_last_exclusive - 1))
            exit()

        self.step = "r"
        orbits = self.success_qa.copy()
        processors = Config.TROPOMI_Processors.copy()
        processors.sort(reverse=True)
        if not(logger is None): logger.info("Start: {} End: {}".format(min(orbits), max(orbits)))

        ncfiles = self.find_highest_configured_procesor(processors)

        for i in range(self.orbit_first_inclusive, self.orbit_last_exclusive):
            if not (self.orbitprocessed_r is None) and i <= self.orbitprocessed_r: continue

            self.orbitprocessed_r = i
            if not i in orbits: 
                logger.info("OMITTED: NC orbit: {} is not included in the initial set".format(i))
                with open(self.pickle, mode="wb") as f:
                    dump(self, f)
                    f.close()
                continue
        
            if not (i in ncfiles): 
                logger.info("OMITTED: NC file for orbit: {} is not available".format(i))
            else:
                f = ncfiles[i]
                processor = TROPOMI.get_processor_version(f)
                filename = path.join(Config.TROPOMI_folder, processor, f)
                tropomi = TROPOMI(filename)
            
                r_CH4_albedo_SWIR, r_CH4_albedo_NIR, r_CH4_areosol_SWIR, r_CH4_areosol_NIR = self.calculate_correlation(tropomi)
                
                # Change NIR is relevant to middle troposphere as it does not resolve differences where the surface temperature 
                # and methane temperature are close. Therefore use only filtering fo SWIR which resolves low level
                # if (0.5 < r_CH4_albedo_SWIR) or (0.5 < r_CH4_albedo_NIR) or (0.5 < r_CH4_areosol_SWIR) or (0.5 < r_CH4_areosol_NIR):
                if (0.5 < r_CH4_albedo_SWIR) or (0.5 < r_CH4_areosol_SWIR):
                    m = "FAIL orbit: {}, R albedo SWIR: {}, R albedo NIR: {}, R aerosol SWIR: {}, R aerosol NIR: {}".format(i,r_CH4_albedo_SWIR, r_CH4_albedo_NIR, r_CH4_areosol_SWIR, r_CH4_areosol_NIR)
                    if (not logger is None): logger.info(m)
                else:
                    m = "PASS orbit: {}, R albedo SWIR: {}, R albedo NIR: {}, R aerosol SWIR: {}, R aerosol NIR: {}".format(i,r_CH4_albedo_SWIR, r_CH4_albedo_NIR, r_CH4_areosol_SWIR, r_CH4_areosol_NIR)
                    if (not logger is None): logger.info(m)
                    self.success_r.append(i)
            
            self.result = self.success_r.copy()
            # Pickle results
            with open(self.pickle, mode="wb") as f:
                dump(self, f)
                f.close()

    def find_highest_configured_procesor(self, processors: list[str]) -> dict[int, str]:
        '''
        Within local directories find the highest configured processor for each orbit.
        At exit dictionary ncfiles contain a mapping orbit to filename.
        
        Parameters
        ----------
        processors: Array of processor versions e.g. ["010301", "010302"]

        Return
        ------
        ncfiles: Dictionary orbit to filename
        '''
        ncfiles = dict()
        for p in processors:
            files = listdir(path.join(Config.TROPOMI_folder, p))
            for f in files:
                if (f.endswith(".nc") and ("L2__CH4____" == TROPOMI.get_productidentifier(f))):
                    orbit = TROPOMI.get_orbit(f)
                    iorbit = int(orbit)
                    if iorbit in ncfiles:
                        pold = TROPOMI.get_processor_version(ncfiles[iorbit])
                        if pold < p:
                            ncfiles[iorbit] = f    
                    else: 
                        ncfiles[iorbit] = f
        return ncfiles

    def calculate_correlation(self, tropomi):
        '''
        Take TROPOMI file and calculate 4 correlation coefficients.
        This method is public to allow re-use by Paper

        Parameters
        ----------
        tropomi - TROPOMI object

        Returns
        -------
        tuple(r_CH4_albedo_SWIR, r_CH4_albedo_NIR, r_CH4_areosol_SWIR, r_CH4_areosol_NIR)
            Four correlation coefficients
        '''
        v_ch4 = tropomi.get_variable_data("methane_mixing_ratio_bias_corrected")[0, :, :]
        v_surface_albedo_SWIR = tropomi.get_variable_data("surface_albedo_SWIR")[0, :, :]
        v_surface_albedo_NIR = tropomi.get_variable_data("surface_albedo_NIR")[0, :, :]
        v_aeorosol_SWIR = tropomi.get_variable_data("aerosol_optical_thickness_SWIR")[0, :, :]
        v_aeorosol_NIR = tropomi.get_variable_data("aerosol_optical_thickness_NIR")[0, :, :]

        l_ch4_albedo_SWIR = []
        l_ch4_albedo_NIR = []
        l_ch4_areosol_SWIR = []
        l_ch4_areosol_NIR = []

        albedo_SWIR = []
        albedo_NIR = []
        areosol_SWIR = []
        areosol_NIR = []
        lt_box_arr = tropomi.get_variable_data("latitude_bounds")[0, :, :, :]
        ln_box_arr = tropomi.get_variable_data("longitude_bounds")[0, :, :, :]
        ln_center_arr = tropomi.get_variable_data("longitude")[0, :, :]
        lt_center_arr = tropomi.get_variable_data("latitude")[0, :, :]
        for s in range(tropomi.scans):
            for p in range(tropomi.pixels):
                if self.config.Algorithm_CSF_Filter_R_domain:
                    if geometry == Geometry.CENTER:
                        ln = ln_center_arr[s,p]
                        lt = lt_center_arr[s,p]
                        p = Point(ln, lt)
                        if not self.config.Algorithm_CSF_domain.contains(p): continue
                    elif geometry == Geometry.CONTAINS:
                        ul = (ln_box_arr[s,p,3], lt_box_arr[s,p,3])
                        bl = (ln_box_arr[s,p,0], lt_box_arr[s,p,0])
                        br = (ln_box_arr[s,p,1], lt_box_arr[s,p,1])
                        ur = (ln_box_arr[s,p,2], lt_box_arr[s,p,2])
                        pixelbox = Polygon([ul,bl,br,ur])
                        if not self.config.Algorithm_CSF_domain.contains(pixelbox): continue
                    elif geometry == Geometry.INTERSECTS:
                        ul = (ln_box_arr[s,p,3], lt_box_arr[s,p,3])
                        bl = (ln_box_arr[s,p,0], lt_box_arr[s,p,0])
                        br = (ln_box_arr[s,p,1], lt_box_arr[s,p,1])
                        ur = (ln_box_arr[s,p,2], lt_box_arr[s,p,2])
                        pixelbox = Polygon([ul,bl,br,ur])
                        if not self.config.Algorithm_CSF_domain.intersects(pixelbox): continue
                    
                c = v_ch4[s,p]                 
                if not(v_ch4.mask[s,p]):
                    if not(v_surface_albedo_SWIR.mask[s,p]):
                        l_ch4_albedo_SWIR.append(c)
                        albedo_SWIR.append(v_surface_albedo_SWIR[s,p])
                            
                    if not(v_surface_albedo_NIR.mask[s,p]):
                        l_ch4_albedo_NIR.append(c)
                        albedo_NIR.append(v_surface_albedo_NIR[s,p])

                    if not(v_aeorosol_SWIR.mask[s,p]):
                        l_ch4_areosol_SWIR.append(c)
                        areosol_SWIR.append(v_aeorosol_SWIR[s,p])
                            
                    if not(v_aeorosol_NIR.mask[s,p]):
                        l_ch4_areosol_NIR.append(c)
                        areosol_NIR.append(v_aeorosol_NIR[s,p])

        r_CH4_albedo_SWIR = TROPOMI_Filter.correlation(l_ch4_albedo_SWIR, albedo_SWIR)
        r_CH4_albedo_NIR = TROPOMI_Filter.correlation(l_ch4_albedo_NIR, albedo_NIR)
        r_CH4_areosol_SWIR = TROPOMI_Filter.correlation(l_ch4_areosol_SWIR, areosol_SWIR)
        r_CH4_areosol_NIR = TROPOMI_Filter.correlation(l_ch4_areosol_NIR, areosol_NIR)
        return r_CH4_albedo_SWIR,r_CH4_albedo_NIR,r_CH4_areosol_SWIR,r_CH4_areosol_NIR

    @staticmethod
    def correlation(l1, l2):
        mean1 = ma.mean(l1)
        mean2 = ma.mean(l2)
        sxy = 0
        sx = 0
        sy = 0
        ln = len(l1)
        for i in range(ln):
            sxy += (l1[i] - mean1) * (l2[i] - mean2)
            sx += (l1[i] - mean1) * (l1[i] - mean1)
            sy += (l2[i] - mean2) * (l2[i] - mean2)
        return fabs(sxy / sqrt(sx) / sqrt(sy))
    
    def list_Properties(self, args):
        '''
        List current status of the filter
        '''
        if args.step == "all" or args.step == "contains":
            print("Filter contains")
            print("\tFirst date inclusive: {} Last date exclusive: {}".format(self.date_first_inclusive, self.date_last_exclusive))
            print("\tDomain bounds: [{0[0]}, {0[1]}, {0[2]}, {0[3]}] ".format(self.config.Algorithm_CSF_domain.bounds))
            print("\tLast date processed: {:04d}-{:02d}-{:02d}".format(self.dateprocessed.year, self.dateprocessed.month, self.dateprocessed.day))
            print("\tFiltered count: {}".format(len(self.success_contains_s3_key)))
            print("\tFiltered orbits:")
            print(self.success_contains)
            print("Filter completed: {}".format(self.completed("contains")))

        if args.step == "all": print("")

        if args.step == "all" or args.step == "qa":
            print("Filter qa")
            print("\tFirst orbit inclusive: {} Last orbit exclusive: {}".format(self.orbit_first_inclusive, self.orbit_last_exclusive))
            print("\tLast orbit processed: {}".format(self.orbitprocessed_qa))
            print("\tFiltered count: {}".format(len(self.success_qa)))
            print("\tFiltered orbits:")
            print(self.success_qa)
            print("Filter completed: {}".format(self.completed("qa")))
        if args.step == "all": print("")

        if args.step == "all" or args.step == "r":
            print("Filter r")
            print("\tFirst orbit inclusive: {} Last orbit exclusive: {}".format(self.orbit_first_inclusive, self.orbit_last_exclusive))
            print("\tLast orbit processed: {}".format(self.orbitprocessed_r))
            print("\tFiltered count: {}".format(len(self.success_r)))
            print("\tFiltered orbits:")
            print(self.success_r)
            print("Filter completed: {}".format(self.completed("r")))

    def orbit_to_date(self, orbit: int) -> datetime:
        '''
        Return endtime of the orbit that is in contains filter
        '''
        l1 = 105
        l2 = 120
        for i in range(len(self.success_contains)):
            if orbit == self.success_contains[i]:
                s = self.success_contains_s3_key[i]
                t = s[l1:l2]
                return datetime(int(t[0:4]), int(t[4:6]), int(t[6:8]), tzinfo = timezone.utc)
        return None

def _handler_clean(args: Namespace):
    '''
    Handler for cleaning up filter
    '''
    Config.create_log("TROPOMI_Filter_clean.log")
    if args.step == "contains" and path.exists(TROPOMI_Filter.pickle):
        remove(TROPOMI_Filter.pickle)
        return

    if args.step == "qa" and path.exists(TROPOMI_Filter.pickle):

        with open(TROPOMI_Filter.pickle, "rb") as f:
            filter: TROPOMI_Filter = load(f)
            f.close()

        filter.orbitprocessed_r = None
        filter.success_r = []
        filter.orbitprocessed_qa = None
        filter.success_qa = []

        with open(TROPOMI_Filter.pickle, "wb") as f:
            dump(filter, f)
            f.close()

        return

    if args.step == "r" and path.exists(TROPOMI_Filter.pickle):

        with open(TROPOMI_Filter.pickle, "rb") as f:
            filter: TROPOMI_Filter = load(f)
            f.close()

        filter.orbitprocessed_r = None
        filter.success_r = []

        with open(TROPOMI_Filter.pickle, "wb") as f:
            dump(filter, f)
            f.close()

        return

def _handler_find(args: Namespace):
    Config.create_log("TROPOMI_Filter_find.log")

    if not path.exists(TROPOMI_Filter.pickle):
        m = "Filter output does not exist: {}".format(TROPOMI_Filter.pickle)
        print(m)
        if not logger is None: logger.error(m)
        return
        
    with open(TROPOMI_Filter.pickle, "rb") as f:
        filter: TROPOMI_Filter = load(f)
        f.close()

    for m in filter.success_contains_s3_key:
        if args.substring in m:
            print(m)

def _handler_list(args: Namespace):
    '''
    Handler for listing filtered orbits 
    '''
    config = Config()

    if not path.exists(TROPOMI_Filter.pickle):
        m = "Filter output does not exist: {}".format(TROPOMI_Filter.pickle)
        print(m)
        return
        
    with open(TROPOMI_Filter.pickle, "rb") as f:
        filter: TROPOMI_Filter = load(f)
        f.close()

    filter.list_Properties(args)

def _handler_run_contains(args: Namespace):
    '''
    Handler for executing this filter
    '''
    Config.create_log("TROPOMI_Filter_run.log")

    # Continue partial calculations
    if path.exists(TROPOMI_Filter.pickle):       
        with open(TROPOMI_Filter.pickle, "rb") as f:
            filter: TROPOMI_Filter = load(f)
            f.close()
        
        if filter.completed("contains"):
            print("Filter contains completed")
            return
        else:
            filter._filter_1_contains(args.username, args.password)
            print("Filter contains completed")
            return
    
    # Start calculations from scratch
    else:
        startdate = datetime(2018,4,30, tzinfo=timezone.utc)
        enddate = datetime(2020,1,1, tzinfo=timezone.utc)
        config = Config()
        if not logger is None: logger.info(config)
        filter = TROPOMI_Filter("contains", config, startdate, enddate)
        filter._filter_1_contains(args.username, args.password)

def _handler_run_qa(args: Namespace):
    '''
    Handler for executing this filter
    '''
    Config.create_log("TROPOMI_Filter_run.log")

    # Continue partial calculations
    if path.exists(TROPOMI_Filter.pickle):       
        with open(TROPOMI_Filter.pickle, "rb") as f:
            filter: TROPOMI_Filter = load(f)
            f.close()

        if not filter.completed("qa"): filter._filter_2_quality()
        print("Filter qa completed")
        return
    
    else:
        m = "Error: You need to run 'python .\\TROPOMI_Filter.py run contains username password' first"
        if not logger is None: logger.error(m)
        print(m)
        return

def _handler_run_r(args: Namespace):
    '''
    Handler for executing this filter
    '''
    Config.create_log("TROPOMI_Filter_run.log")

    # Continue partial calculations
    if path.exists(TROPOMI_Filter.pickle):       
        with open(TROPOMI_Filter.pickle, "rb") as f:
            filter: TROPOMI_Filter = load(f)
            f.close()

        if not filter.completed("r"): filter._filter_3_correlation()
        print("Filter r completed")
        return
    else:
        m = "Error: You need to run 'python .\\TROPOMI_Filter.py run contains username password' first. Then you need to run 'python .\\TROPOMI_Filter.py run qa'"
        if not logger is None: logger.error(m)
        print(m)
        return

if __name__ == '__main__': 

    help_source = "Source of emissions"
    parser = ArgumentParser(prog="TROPOMI_Filter")
    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_clean = subparsers.add_parser("clean", help="Clean filter status.")
    parser_clean.add_argument("step", choices=["contains", "qa", "r"], help="contains - removes r, qa and contains results; qa - removes r and qa; r - removes r")
    parser_clean.set_defaults(func=_handler_clean)

    parser_find = subparsers.add_parser("find", help="Find filtered file names containing susbtring.")
    parser_find.add_argument("substring", type=str, help="Can be a date as YYYYMMDD or orbit as 5 digit string")
    parser_find.set_defaults(func=_handler_find)

    parser_list = subparsers.add_parser("list", help="List filtered files containing the source domain.")
    parser_list.add_argument("step", choices=["all", "contains", "qa", "r"], help="filter step")
    parser_list.set_defaults(func=_handler_list)

    parser_run = subparsers.add_parser("run", help="Run filters using THREDDS files for the source domain.")
    subparsers_run = parser_run.add_subparsers(help="subcommand help", required=True)

    parser_run_contains = subparsers_run.add_parser("contains", help="Filter files containing domain within specified range of orbits")
    parser_run_contains.add_argument("username", type=str, help="Copernicus EU username")
    parser_run_contains.add_argument("password", type=str, help="Copernicus EU password")
    parser_run_contains.set_defaults(func=_handler_run_contains)

    parser_run_qa = subparsers_run.add_parser("qa", help="Filter files containing more than 500 qa = 1")
    parser_run_qa.set_defaults(func=_handler_run_qa)

    parser_run_r = subparsers_run.add_parser("r", help="Filter files proocessed by qa filter. Based on correlation to albedo and aeorosol.")
    parser_run_r.set_defaults(func=_handler_run_r)


    args = parser.parse_args()
    
    args.func(args)