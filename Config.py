from Geometry import Geometry
from logging import basicConfig, INFO
from Mask import Mask
from os import mkdir, path
from shapely import box, prepare
class Config: 

    #### Start Configure Local Directories
    '''
    Local directory containing TROPOMI files.
    '''        
    Root_folder = path.join("..", "Test")
    BoM_folder = path.join(Root_folder, "BoM")
    CSF_folder = path.join(Root_folder, "CSF")
    ERA5_folder = path.join(Root_folder, "ERA5")
    HYSPLIT_folder = path.join(Root_folder, "HYSPLIT")
    LOG_folder = path.join("Log")
    TROPOMI_folder = path.join(Root_folder, "TROPOMI")

    # Start "010202", "010300", "010302", "010400", "020200", "020301"
    Paper_folder = path.join(Root_folder, "Paper", "Original")
    TROPOMI_Processors = ["010202", "010300", "010301", "010302"] 
    # End "010202", "010300", "010302", "010400", "020200", "020301"

    # Start "020400"

    # Paper_folder = path.join("Root_folder, "Paper", "020400")
    # TROPOMI_Processors = ["020400"] # Files post processed in November 2022

    # End "020400"

    #### End Configure Local Directories
        
    def __init__(self):
        '''
        Configuration of software.
        1. Configure Algorithm CSF
        2. Configure TROPOMI
        3. Configure local directories storing input and output
        '''
        #### Start Section Algorithm CSF replicating Sadaverte 2021
        
        '''
        Geometry used for domain median and upwind box average background calculations
        '''
        self.Algorithm_CSF_background_geometry = Geometry.CONTAINS
        
        '''
        Minimum number of pixels in the upwind box to use the average following Sadavarte et all 2021
        '''
        self.Algorithm_CSF_background_minimum_count = 20

        '''
        Domain configuration following Sadavarte et all 2021
        '''
        self.Algorithm_CSF_domain_xmax = 150.0
        self.Algorithm_CSF_domain_xmin = 146.0
        self.Algorithm_CSF_domain_ymax = -20.0
        self.Algorithm_CSF_domain_ymin = -24.0

        self.Algorithm_CSF_domain = box(
            self.Algorithm_CSF_domain_xmin, 
            self.Algorithm_CSF_domain_ymin, 
            self.Algorithm_CSF_domain_xmax, 
            self.Algorithm_CSF_domain_ymax
        )

        '''
        Geometry used for downwind box shape calculations
        '''
        self.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS

        '''
        Mask used for downwind box shape calculations
        '''
        self.Algorithm_CSF_downwindbox_mask = Mask.NONE

        '''
        Minimum enhancement increment in ppb
        '''
        self.Algorithm_CSF_downwindbox_enhancement_delta = 5

        '''
        The first orbit available om NCI THREDDS server 30/Apr/2018. Used to narrow filtering
        '''
        self.Algorithm_CSF_first_orbit_inclusive = "02818"

        # The first orbit on January first 2020 available at NCI THREDDS server. Used to narrow filtering
        self.Algorithm_CSF_last_orbit_exclusive = "11487"

        # Background calculation method used for calculations of transect enhancement
        # "algorithm" - Use background value as determined by background Sadavarte algorithm
        # "domain" - Force to use domain median ignoring algorithm. Note that algorithm is still used for downwind box shape determination
        # "manual" - Use value manually specfied at the command line
        # "minimum" - Use the minimum value in the downwind box.
        # "upwind" - Force to use the average of the upwind box
        self.Algorithm_CSF_transect_background = "algorithm"

        # Value of background used for manual option
        self.Algorithm_CSF_transect_background_value: float = None

        # For calculation of average emissions filter out orbits where transect intersect fewer then 15 positive pixels 
        self.Algorithm_CSF_transect_positive_minimum_count: int = 15

        # For calculation of average emissions filter out orbits where transect intersect more then 70 valid pixels
        # This is likely overesitmation of plume
        self.Algorithm_CSF_transect_valid_maximumu_count: int = 70

        # If true calculate correlations based over algorithm domain. If false use correlations over the whole orbit
        self.Algorithm_CSF_Filter_R_domain: bool = False

        # End Section Algorithm CSF

        prepare(self.Algorithm_CSF_domain)

    def __eq__(self, other):
        '''
        Override of the equality for configurations
        '''
        if not isinstance(other, Config): return False
        if not self.Algorithm_CSF_background_geometry == other.Algorithm_CSF_background_geometry: return False
        if not self.Algorithm_CSF_background_minimum_count == other.Algorithm_CSF_background_minimum_count: return False
        if not self.Algorithm_CSF_domain_xmax == other.Algorithm_CSF_domain_xmax: return False
        if not self.Algorithm_CSF_domain_xmin == other.Algorithm_CSF_domain_xmin: return False
        if not self.Algorithm_CSF_domain_ymax == other.Algorithm_CSF_domain_ymax: return False
        if not self.Algorithm_CSF_domain_ymin == other.Algorithm_CSF_domain_ymin: return False
        if not self.Algorithm_CSF_downwindbox_geometry == other.Algorithm_CSF_downwindbox_geometry: return False
        if not self.Algorithm_CSF_downwindbox_mask == other.Algorithm_CSF_downwindbox_mask: return False
        if not self.Algorithm_CSF_downwindbox_enhancement_delta == other.Algorithm_CSF_downwindbox_enhancement_delta: return False
        if not self.Algorithm_CSF_first_orbit_inclusive == other.Algorithm_CSF_first_orbit_inclusive: return False
        if not self.Algorithm_CSF_last_orbit_exclusive == other.Algorithm_CSF_last_orbit_exclusive: return False
        if not self.Algorithm_CSF_transect_background == other.Algorithm_CSF_transect_background: return False
        if not self.Algorithm_CSF_transect_background_value == other.Algorithm_CSF_transect_background_value: return False
        if not self.Algorithm_CSF_Filter_R_domain == other.Algorithm_CSF_Filter_R_domain: return False
        return True

    def __str__(self) -> str:
        '''
        List configuration
        '''
        m = []
        m.append("Config")
        m.append("\tDomain: Lower left corner: [{}, {}] Upper right corner: [{}, {}]".format(
            self.Algorithm_CSF_domain_xmin,
            self.Algorithm_CSF_domain_ymin, 
            self.Algorithm_CSF_domain_xmax, 
            self.Algorithm_CSF_domain_ymax
        ))
        m.append("\tCSF Background geometry: {}".format(self.Algorithm_CSF_background_geometry.name))
        m.append("\tCSF Background minimum count: {}".format(self.Algorithm_CSF_background_minimum_count))
        m.append("\tCSF Downwind box geometry: {}".format(self.Algorithm_CSF_downwindbox_geometry.name))
        m.append("\tCSF Downwind box mask: {}".format(self.Algorithm_CSF_downwindbox_mask.name))
        m.append("\tCSF Downwind box enhancement delta: {} (ppb)".format(self.Algorithm_CSF_downwindbox_enhancement_delta))
        m.append("\tCSF First orbit inclusive: {}".format(self.Algorithm_CSF_first_orbit_inclusive))
        m.append("\tCSF Last orbit exclusive: {}".format(self.Algorithm_CSF_last_orbit_exclusive))
        m.append("\tCSF Background algorithm: {}".format(self.Algorithm_CSF_transect_background))
        if not self.Algorithm_CSF_transect_background_value is None: m.append("\tCSF Background Manual value: {} (ppb)".format(self.Algorithm_CSF_transect_background_value))
        m.append("\tTROPOMI processors: {}".format(",".join(Config.TROPOMI_Processors)))
        return "\n".join(m)
        
    @staticmethod
    def _create_directory_root():
        p_root = path.join(Config.Root_folder)
        if not (path.exists(p_root)) :  mkdir(p_root, 0o777)

    @staticmethod
    def create_directory_structure_BOM(case_name: str, orbit: str):
        ''' 
        Set up directory structure for case study

        Parameters
        ----------
        casename: str
            Name of folder for this study case
        orbit: str
            Orbit identifier
        '''
        Config._create_directory_root()
        p_bom = path.join(Config.BoM_folder)
        if not (path.exists(p_bom)) :  mkdir(p_bom, 0o777)

        if not(case_name is None) and 0 < len(case_name):
            p_bom_casename = path.join(p_bom, case_name)
            if not (path.exists(p_bom_casename)) :  mkdir(p_bom_casename, 0o777)
            
            if not (orbit is None) and 0 < len(orbit):
                p_data_casename_orbit = path.join(p_bom_casename, orbit)
                if not (path.exists(p_data_casename_orbit)) :  mkdir(p_data_casename_orbit, 0o777)
        return None

    @staticmethod
    def create_directory_structure_CSF(case_name: str):
        ''' 
        Set up directory structure for case study

        Parameters
        ----------
        casename: str
            Name of folder for this study case
        '''
        Config._create_directory_root()
        p_csf = path.join(Config.CSF_folder)

        if not (path.exists(p_csf)) :  mkdir(p_csf, 0o777)

        if not(case_name is None) and 0 < len(case_name):
            p_csf_casename = path.join(p_csf, case_name)
            if not (path.exists(p_csf_casename)) :  mkdir(p_csf_casename, 0o777)                
        return None

    @staticmethod
    def create_directory_structure_ERA5(case_name: str):
        ''' 
        Set up directory structure for case study

        Parameters
        ----------
        casename: str
            Name of folder for this study case
        '''
        Config._create_directory_root()
        p_era5 = path.join(Config.ERA5_folder)
        if not (path.exists(p_era5)) :  mkdir(p_era5, 0o777)

        if not(case_name is None) and 0 < len(case_name):
            p_era5_casename = path.join(p_era5, case_name)
            if not (path.exists(p_era5_casename)) :  mkdir(p_era5_casename, 0o777)
            
        return None

    @staticmethod
    def create_directory_structure_HYSPLIT(case_name: str):
        ''' 
        Set up directory structure for HYSPLIT files
        '''
        Config._create_directory_root()
        p_hysplit = path.join(Config.HYSPLIT_folder)
        p_hysplit_casename = path.join(Config.HYSPLIT_folder, case_name)

        if not (path.exists(p_hysplit)) :  mkdir(p_hysplit, 0o777)
        if not (path.exists(p_hysplit_casename)) :  mkdir(p_hysplit_casename, 0o777)

        return None

    @staticmethod
    def create_directory_structure_TROPOMI(processor: str):
        ''' 
        Set up directory structure for case study

        Parameters
        ----------
        orbit: str
            Orbit identifier
        processor: str
            Processor identifier
        '''
        Config._create_directory_root()

        p_tropomi = path.join(Config.TROPOMI_folder)
        if not (path.exists(p_tropomi)) :  mkdir(p_tropomi, 0o777)
            
        if not (processor is None) and 0 < len(processor):
            p_tropomi_processor = path.join(p_tropomi, processor)
            if not (path.exists(p_tropomi_processor)) :  mkdir(p_tropomi_processor, 0o777)
                
        return None

    @staticmethod
    def create_log(filename: str):
        p_Log = path.join(Config.LOG_folder)
        if not (path.exists(p_Log)) :  mkdir(p_Log, 0o777)
        FORMAT = '%(asctime)s %(module)s %(funcName)s %(message)s'
        path_LOGGER = path.join(p_Log, filename)
        basicConfig(filename=path_LOGGER, level=INFO, format=FORMAT)
        return None
 