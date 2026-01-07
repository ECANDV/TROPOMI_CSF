from Algorithm_CSF import Algorithm_CSF, EnhancementLength, EnhancementWidth, Geometry, Mask, Rotation, Status, Transect
from Algorithm_HYSPLIT import Algorithm_HYSPLIT, Algorithm_HYSPLIT_from_trajectory
from Algorithm_IME import Algorithm_IME, Algorithm_IME_factory
from Algorithm_TM import Algorithm_TM, Algorithm_TM_factory
from argparse import ArgumentParser
from BOM_AWS import create_aws, BOM_AWS_from_csv, path_Moranbah_AWS_CSV
from BOM_MSLP import chart_MSLP
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from Config import Config
from Constants import GRAVITY
from datetime import datetime, timedelta, timezone
from ERA5 import ERA5, ERA5_SingleLevel, ERA5_PressureLevels
from Geometry import Geometry
from logging import getLogger
from Mask import Mask
from math import cos, pi, sin, sqrt
from Meteorology import Meteorology
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from netCDF4 import Dataset
from numpy import average, linspace, ma, max, min, std
from os import listdir, path, remove
from pickle import load
from shapely import box, Point, Polygon
from Source import create_source, Source
from SRTM import path_9s_geotiff, SRTM_from_geotiff
from trajectory import drawTrajectory
from TROPOMI import TROPOMI, TROPOMI_for_orbit
from TROPOMI_Filter import TROPOMI_Filter

logger = getLogger(__name__)

# Source: Methane Emissions from Super-emitting Coal Mines in Australia quantified using TROPOMI Satellite Observations
# Emissions by Pankaj Sadavarte*,†, Sudhanshu Pandey†, Joannes D. Maasakkers†, Alba Lorente†, 
# Tobias Borsdorff†, Hugo Denier van der Gon‡, Sander Houweling†,§, Ilse Aben†
csv_hailcreek = ["Details,Source 1"]
csv_hailcreek.append("location, Hail Creek")
csv_hailcreek.append("mine type, surface")
csv_hailcreek.append("mining method, Dragline; truck and shovel")
csv_hailcreek.append("total raw coal production in million tonnes 2018 - 19, 7.7")
csv_hailcreek.append("total raw coal production in million tonnes 2019 - 20, 5.8")
csv_hailcreek.append("longitude, 148.380°E")
csv_hailcreek.append("latitude, 21.490°S")
csv_hailcreek.append("number of clear-sky observations in TROPOMI, 32")
csv_hailcreek.append("annual emissions using the CSF method (Gg/year) [mu ± 2sigma], 230 ± 50")

# Source: Methane Emissions from Super-emitting Coal Mines in Australia quantified using TROPOMI Satellite Observations
# Emissions by Pankaj Sadavarte*,†, Sudhanshu Pandey†, Joannes D. Maasakkers†, Alba Lorente†, 
# Tobias Borsdorff†, Hugo Denier van der Gon‡, Sander Houweling†,§, Ilse Aben†
csv_ef = ["Source,Reference,Emission year,IEF (g/kg),Comment"]
csv_ef.append("Australia,EDGAR v4.3.2 Table 1.B.1 CRF NIR 2020,2012,0.55,Split EDGAR underground and surface based on 2012 ratio NIR 2017.EDGARv4.3.2 CH4 emissions for 2012 and used national raw coalproduction from common reporting format table for 2012.")
csv_ef.append("Australia,Table 1.B.1 CRF NIR 2020,2012,0.45,Table 1.B.1")
csv_ef.append("Australia,Table 1.B.1 CRF NIR 2020,2017,0.47,Table 1.B.1")
csv_ef.append("Australia,Table 1.B.1 CRF NIR 2020,2018,0.53,Table 1.B.1")
csv_ef.append("Queensland,NIR (2020),2018,0.73,Using emissions and raw coal production details from commonreporting format table for respective year.");
csv_ef.append("Hail Creek, Sadavarte et all. 2021,2018,0.73 g_kg-1,Reconstructed bottom-up")
csv_ef.append("Global, Kholod et all. 2020,,2.03 - 3.38,below 200m")
csv_ef.append("Hail Creek, Sadavarte et all. 2021,34.12, CSF")

# Source IPCC 2006 Guidelines Volume 2 Chapter 4 page 4.18 https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/2_Volume2/V2_4_Ch4_Fugitive_Emissions.pdf
csv_ef_ipcc = ["Source,Reference,IEF (g/kg),Comment"]
csv_ef_ipcc.append("Global,p.4.18 (IPCC 2006),0.20,Low Emission Factor. Value in (Sadavarte, et al., Methane Emissions from Super-emitting Coal Mines in Australia quantified using TROPOMI Satellite Observations. Supporting Information, 2021) marked as 0 - 200m.")
csv_ef_ipcc.append("Global,p.4.18 (IPCC 2006),0.88,Average Emission Factor. Value in (Sadavarte, et al., Methane Emissions from Super-emitting Coal Mines in Australia quantified using TROPOMI Satellite Observations. Supporting Information, 2021) is slightly higher than conversion and marked as 200-400m")
csv_ef_ipcc.append("Global,p.4.18 (IPCC 2006),1.49,High Emission Factor. Value in (Sadavarte, et al., Methane Emissions from Super-emitting Coal Mines in Australia quantified using TROPOMI Satellite Observations. Supporting Information, 2021) is slightly higher than conversion and marked as 400+ m")
csv_ef.append("Source: p.4.18 (IPCC 2006)")

# Source Qld Gov Open Data https://www.data.qld.gov.au/dataset/annual-coal-statistics/resource/d22a8d8b-7c00-42d2-884a-c438d51cefc3
csv_activity_cy = ["Calendar Year,Total Net Output (tonnes),Total Gross Raw Output (tonnes)"]
csv_activity_cy.append("2022,,9833100")
csv_activity_cy.append("2021,,11629180")
csv_activity_cy.append("2020,,9849473")
csv_activity_cy.append("2019,8751078,,")
csv_activity_cy.append("2018,9028926,,")
csv_activity_cy.append("2017,9500316,,")

# Source https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
csv_activity_fy = ["Financial year,Coal type, Activity Discards,Gross Raw Feed,Gross Raw Output,Net Output"]
csv_activity_fy.append("FY2025,Coal,,16154750,10807512,")
csv_activity_fy.append("FY2025,Coking,,,,4452342")
csv_activity_fy.append("FY2025,Thermal,,,,3469428")
csv_activity_fy.append("FY2025,Wash Discard,8803589,,,")
csv_activity_fy.append("FY2024,Coal,,16575361,9592586,")
csv_activity_fy.append("FY2024,Coking,,,,3895983")
csv_activity_fy.append("FY2024,Thermal,,,,3977955")
csv_activity_fy.append("FY2024,Wash Discard,9598043,,,")
csv_activity_fy.append("FY2023,Coal,,16835413,10092999,")
csv_activity_fy.append("FY2023,Coking,,,,4545708")
csv_activity_fy.append("FY2023,Thermal,,,,3082427")
csv_activity_fy.append("FY2023,Wash Discard,9389503,,,")
csv_activity_fy.append("FY2022,Coal,,17436110,10783931,")
csv_activity_fy.append("FY2022,Coking,,,,5140618")
csv_activity_fy.append("FY2022,Thermal,,,,2658374")
csv_activity_fy.append("FY2022,Wash Discard,9766081,,,")
csv_activity_fy.append("FY2021,Coal,,13024250,10594483,")
csv_activity_fy.append("FY2021,Coking,,,,4698426")
csv_activity_fy.append("FY2021,Thermal,,,,1024118")
csv_activity_fy.append("FY2021,Wash Discard,7223108,,,")
csv_activity_fy.append("FY2020,Coal,,14878458,9468696,")
csv_activity_fy.append("FY2020,Coking,,,,5441586")
csv_activity_fy.append("FY2020,Thermal,,,,3014628")
csv_activity_fy.append("FY2020,Wash Discard,3819050,,,")
csv_activity_fy.append("FY2019,Coal,,15616001,7660496,")
csv_activity_fy.append("FY2019,Coking,,,,5535276")
csv_activity_fy.append("FY2019,Thermal,,,,3612680")
csv_activity_fy.append("FY2019,Wash Discard,2386331,,,")
csv_activity_fy.append("FY2018,Coal,,10000076,10203613,")
csv_activity_fy.append("FY2018,Coking,,,,5356850")
csv_activity_fy.append("FY2018,Thermal,,,,4175818")
csv_activity_fy.append("FY2018,Wash Discard,967078,,,")
csv_activity_fy.append("FY2017,Coal,,9166559,9195752,")
csv_activity_fy.append("FY2017,Coking,,,,5178543")
csv_activity_fy.append("FY2017,Thermal,,,,4032202")
csv_activity_fy.append("FY2017,Wash Discard,1303906,,,")
csv_activity_fy.append("FY2016,Coal,,10439914,10228814,")
csv_activity_fy.append("FY2016,Coking,,,,6071692")
csv_activity_fy.append("FY2016,Thermal,,,,3455348")
csv_activity_fy.append("FY2016,Wash Discard,1523574,,,")
csv_activity_fy.append("FY2015,Coal,,11709624,11614548,")
csv_activity_fy.append("FY2015,Coking,,,,6392465")
csv_activity_fy.append("FY2015,Thermal,,,,3110521")
csv_activity_fy.append("FY2015,Wash Discard,2287804,,,")
csv_activity_fy.append("FY2014,Coal,,11726258,11911568,")
csv_activity_fy.append("FY2014,Coking,,,,6748102")
csv_activity_fy.append("FY2014,Thermal,,,,989488")
csv_activity_fy.append("FY2014,Wash Discard,3988668,,,")
csv_activity_fy.append("FY2013,Coal,,11804946,12063330,")
csv_activity_fy.append("FY2013,Coking,,,,7068582")
csv_activity_fy.append("FY2013,Wash Discard,4736364,,,")
csv_activity_fy.append("FY2012,Coal,,12767641,12817071,")
csv_activity_fy.append("FY2012,Coking,,,,7485638")
csv_activity_fy.append("FY2012,Wash Discard,5282003,,,")
csv_activity_fy.append("FY2011,Coal,,12025272,11983327,")
csv_activity_fy.append("FY2011,Coking,,,,6878982")
csv_activity_fy.append("FY2011,Wash Discard,5146290,,,")
csv_activity_fy.append("FY2010,Coal,,12349908,12554146,")
csv_activity_fy.append("FY2010,Coking,,,,6631738")
csv_activity_fy.append("FY2010,Wash Discard,5718170,,,")

# Source: Compilation of data from sources listed in the first column
csv_emissions = ["Source,CH4 (Gg/year),Period,Method"]
csv_emissions.append("Sadavarte et all. 2021,230.0,Apr 2018 - Dec 2019,TROPOMI -> CSF")
csv_emissions.append("Palmer et all. 2021,42.9,Apr 2018 - Dec 2018,TROPOMI -> GEOS-Chem -> Ensemble Kalman Filter")
csv_emissions.append("Open Methane (beta),19.5,Jan 2023 - Jun 2023,TROPOMI -> WRF -> CMAQ")
csv_emissions.append("Borchardt et all. 2025,219,Apr 2018 - Dec 2023,TROPOMI -> WFMD -> CSF ")
csv_emissions.append("Borchardt et all. 2025,122.6,31 May 2022,Aircraft")
csv_emissions.append("Borchardt et all. 2025,122.6,3 Jun 2022,Aircraft")
csv_emissions.append("Borchardt et all. 2025,84,28 Sep 2023	Aircraft")
csv_emissions.append("Borchardt et all. 2025,99,29 Sep 2023,Aircraft")

def __utility_file_does_not_exist(f: str) -> None:
    '''
    Handler for missing file
    '''        
    m = "Error: Required file does not exists: {}".format(f)
    print(m)
    if not logger is None: logger.error(m)
    return

def __utility_rotate(x, y, c_x, c_y, theta):
    '''
    Anticlockwise rotate a point (x,y) around centre (c_x, c_y) by angle theta in radians
    '''
    return (
        (x - c_x) * cos(theta) - (y - c_y) * sin(theta) + c_x, 
        (x - c_x) * sin(theta) + (y - c_y) * cos(theta) + c_y
    )

def __utility_table_output_csv(fileout, csv):
    '''
    Local output for table 
    '''
    if not fileout is None:
        with open(fileout, "wt") as f:
            f.write("\n".join(csv))
            f.close()
        print("Table generated: {}".format(fileout))
    else:
        for m in csv: print(m)

def _chart_AWS_windrose(fileout:str) -> None:
    '''
    Chart 
    '''
    if not path.exists(path_Moranbah_AWS_CSV) :
        __utility_file_does_not_exist(path_Moranbah_AWS_CSV)
        return

    aws = BOM_AWS_from_csv("Moranbah AWS", path_Moranbah_AWS_CSV)
    aws.plot_wind_rose(fileout)
    return

def _chart_HailCreek_Activity(fileout:str) -> None:
    ds = csv_activity_fy[1:]
    l = len(ds)
    coking_net_dict = dict()
    gross_feed_dict = dict()
    gross_raw_dict = dict()
    thermal_net_dict = dict()
    indexes=[]
    for i in range(l):
        arr = ds[i].split(",")
        if not arr[0] in indexes:
            indexes.append(arr[0])
        
        if arr[1] == "Coking":            
            coking_net_dict[arr[0]] = arr[5]

        if arr[1] == "Thermal":
            thermal_net_dict[arr[0]] = arr[5]

        if arr[1] == "Coal":
            gross_feed_dict[arr[0]] = arr[3]
            gross_raw_dict[arr[0]] = arr[4]

    indexes.sort()

    x = []
    coking_net_list = []
    gross_feed_list = []
    gross_raw_list = []
    net_list = []
    thermal_net_list = []
    labels = []
    i = 0    
    for lab in indexes:
        x.append(i)
        labels.append(i + 2010)
        coking_net_list.append(float(coking_net_dict[lab])/1000000 if lab in coking_net_dict else None)
        gross_feed_list.append(float(gross_feed_dict[lab])/1000000 if lab in gross_feed_dict else None)
        gross_raw_list.append(float(gross_raw_dict[lab])/1000000 if lab in gross_raw_dict else None)
        thermal_net_list.append(float(thermal_net_dict[lab])/1000000 if lab in thermal_net_dict else None)
        net_list.append(coking_net_list[i] + (thermal_net_list[i] if not thermal_net_list[i] is None else 0))
        i += 1

    fig, ax = plt.subplots(figsize=(10,5))
    plt.title("Hail Creek Mine, Raw Feed, Raw Output, Net Output.\nSource: Qld Government Open Data")
    ax.set_xlabel("Financial Year")
    ax.set_xlim(xmin=0, xmax = 15)
    ax.set_xticks(x, labels, minor=False)
    ax.set_ylabel("Output (Mt)")
    ax.set_ylim(ymin = 0, ymax = 18)
    ax.plot(x, net_list, linewidth = 2, color="green", label="Net Output")
    ax.plot(x, gross_feed_list, linewidth = 2, color="red", label="Gross Feed")
    ax.plot(x, gross_raw_list, linewidth = 2, color="black", label="Gross Raw")
    plt.legend()
        
    if not fileout is None:
        plt.savefig(fileout)
        print("Chart generated: {}".format(fileout))
    else:
        plt.show()
    
    plt.close("all")

    return

def _chart_CSF_background_box_combined(orbit: str, processor: str, fileout:str):
    '''
    Chart background box used by algorithm of Sadavarte
    '''
    files = []
    config1 = Config()
    config1.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    
    f1 = Algorithm_CSF.get_picklename(config1, "HailCreek", orbit, processor)
    if not path.exists(f1):
        __utility_file_does_not_exist(f1)
        return
    files.append(f1)

    config2 = Config()
    config2.Algorithm_CSF_background_geometry = Geometry.CENTER
    f2 = Algorithm_CSF.get_picklename(config2, "HailCreek", orbit, processor)
    if not path.exists(f2):
        __utility_file_does_not_exist(f2)
        return
    files.append(f2)

    config3 = Config()
    config3.Algorithm_CSF_background_geometry = Geometry.INTERSECTS
    f3 = Algorithm_CSF.get_picklename(config3, "HailCreek", orbit, processor)
    if not path.exists(f3):
        __utility_file_does_not_exist(f3)
        return
    files.append(f3)

    # Note that background calculation does not depend on enhancement so use all    
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10,5), layout="constrained", sharey=True, subplot_kw=dict(projection=ccrs.PlateCarree()))
    plateCarree = ccrs.PlateCarree()
    geodetic = ccrs.Geodetic()

    scans = []
    lats = []
    lons = []
    maxlat = []
    maxlon = []
    minlat = []
    minlon = []
    v_max = []
    v_min = []

    # Determine common box
    for filename in files:
        if not path.exists(filename):
            message = "File not found: {}".format(filename)
            print("File not found: {}".format(filename))
            if not logger is None: logger.error(message)
            continue

        with open(filename, mode="rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()

        title = "Upwind background box. Orbit: {} Processor: {}".format(sad.tropomi.orbit, sad.tropomi.processor_version)

        [_, _, _, _, _, _, ln, lt, scanCH4] = sad.tropomi.narrow_to_domain(sad.config.Algorithm_CSF_background_geometry, sad.upwind_box)
        scans.append(scanCH4)
        lons.append(ln)
        lats.append(lt)
        maxlat.append(max(sad.upwind_box.exterior.coords.xy[1]))
        minlat.append(min(sad.upwind_box.exterior.coords.xy[1]))
        maxlon.append(max(sad.upwind_box.exterior.coords.xy[0]))
        minlon.append(min(sad.upwind_box.exterior.coords.xy[0]))
        v_min.append(ma.min(scanCH4))
        v_max.append(ma.max(scanCH4))

    group_lat_max = max(maxlat)
    group_lat_min = min(minlat)
    group_lon_max = max(maxlon)
    group_lon_min = min(minlon)
    group_v_max = max(v_max)
    group_v_min = max(v_min)

    for filename in files:
        if not path.exists(filename):
            message = "File not found: {}".format(filename)
            print("File not found: {}".format(filename))
            if not logger is None: logger.error(message)
            continue

        with open(filename, mode="rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()

        n = sad.config.Algorithm_CSF_background_geometry.value - 1
        ax = axs.flat[n]
        ax.set_extent(extents=[group_lon_min - 0.1, group_lon_max + 0.1, group_lat_min - 0.1, group_lat_max + 0.1], crs=plateCarree)
        cs = ax.pcolormesh(lons[n], lats[n], scans[n], shading="nearest", vmin=group_v_min, vmax=group_v_max, cmap=plt.cm.rainbow, transform=ccrs.PlateCarree())

        Source.plot_source_in_extent(sad.tropomi_source_date, group_lon_min, group_lat_min, group_lon_max, group_lat_max, ax)
        ax.set_title(sad.config.Algorithm_CSF_background_geometry.name, fontsize="small", loc="center")

        # Plot point 0.1 degrees upwind from Hail Creek
        drawTrajectory(ax, geodetic, plateCarree, sad.upwind_box.exterior.coords)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.left_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'gray'}
    gl.ylabel_style = {'size': 10, 'color': 'gray'} 

    fig.suptitle(title)
    cbar = fig.colorbar(cs, ax=axs, orientation="horizontal", shrink=0.6, pad=0.075)
    cbar.ax.set_xlabel("CH4 ppb")

    if not fileout is None:
        plt.savefig(fileout)
        print("Chart generated: {}".format(fileout))
    else:
        plt.show()
    plt.close("all")
    return

def _chart_CSF_elements(source:str, orbit: str, processor: str, fileout:str) -> None:
    '''
    Chart all elemnents for current configuration of run for orbit
    '''
    config = Config()    
    filename = Algorithm_CSF.get_picklename(config, source, orbit, processor)
    if not path.exists(filename):
        print("File {} does not exist. Terminating".format(filename))
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()    
    sad.chart_lonlat_all_elements(fileout)
    return

def _chart_CSF_transects(source:str, orbit: str, processor: str, fileout: str) -> None:
    '''
    List transects for current configuration of run for orbit
    '''
    config = Config()    
    filename = Algorithm_CSF.get_picklename(config, source, orbit, processor)
    if not path.exists(filename):
        print("File {} does not exist. Terminating".format(filename))
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()
        sad.chart_transects(fileout)
    
    return

def _chart_CSF_positive(source: str, orbit:str, processor: str, fileout:str):
    config = Config()
    p = Algorithm_CSF.get_picklename(config, source, orbit, processor)
    if not path.exists(p):
        __utility_file_does_not_exist(p)
        return

    with open(p, "rb") as f:
        sad:Algorithm_CSF = load(f)
        sad.chart_lonlat_positive(fileout)

    return

def _chart_CSF_positive_config(source: str, orbit:str, processor: str, config: Config, fileout:str) -> None:
    p = Algorithm_CSF.get_picklename(config, source, orbit, processor)
    if not path.exists(p):
        __utility_file_does_not_exist(p)
        return

    with open(p, "rb") as f:
        sad:Algorithm_CSF = load(f)
        sad.chart_lonlat_positive(fileout)

    return 

def _chart_CSF_valid(source: str, orbit:str, processor: str, fileout: str) -> None:
    config = Config()
    p = Algorithm_CSF.get_picklename(config, source, orbit, processor)
    if not path.exists(p):
        __utility_file_does_not_exist(p)
        return

    with open(p, "rb") as f:
        sad:Algorithm_CSF = load(f)
        sad.chart_lonlat_valid(fileout)

    return

def _chart_HYSPLIT_time_slice(
                         sourcename: str, 
                         orbit: str,
                         processor: str, 
                         model: str, 
                         direction: str, 
                         filter_start_start: datetime, 
                         filter_start_end: datetime, 
                         filter_end_start: datetime, 
                         filter_end_end: datetime, 
                         slice_datetime: datetime,
                         fileout: str) -> None:
    '''
    Chart time slices of filtered trajectories
    '''
    sourcepath = path.join(Config.HYSPLIT_folder, sourcename)
    filenames = listdir(sourcepath)
    pattern = "HYSPLIT_{}_{}_".format(model, direction)
    source = create_source(sourcename)

    ah = Algorithm_HYSPLIT(source, model, slice_datetime)
    for f in filenames:
        if not f.startswith(pattern): continue
        if not logger is None: logger.info("Processing file {}".format(f))
        trajectories = Algorithm_HYSPLIT_from_trajectory(source, model, slice_datetime, path.join(sourcepath, f)).trajectories
        filtered = []
        for t in trajectories:
            s_dt = t.start_datetime()
            e_dt = t.end_datetime()
            if filter_start_start <= s_dt and s_dt <= filter_start_end and filter_end_start <= e_dt and e_dt <= filter_end_end:
                filtered.append(t)
                if not logger is None: logger.info("Adding trajectory from file: {}".format(f))
        ah.addTrajectories(filtered)

    if 0 == len(ah.trajectories):
        if not logger is None: logger.error("No HYSPLIT trajectory files found")
        print("No HYSPLIT trajectory files found")
        return None

    ah.chart_time_slices(orbit, processor, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)

    return None

def _chart_HYSPLIT_HailCreek_SingleTrajectory(sourcename: str, orbit:str, processor:str, model: str, emission_start_date, emission_end_date, chart_datetime, fileout:str) -> None:
    '''
    Draw one forward and one backward trajectory which respectively originates and ends at HailCreek 20190914 23Z
    '''
    source = create_source(sourcename)
    t1 = emission_start_date.strftime("%Y%m%d%H")
    t2 = emission_end_date.strftime("%Y%m%d%H")
    t3 = chart_datetime.strftime("%Y%m%d%H")

    ah = Algorithm_HYSPLIT(source, model, chart_datetime)
    t_backward = path.join(Config.HYSPLIT_folder, sourcename, "HYSPLIT_GFSQ_B_2019091404_2019091423.txt")
    t_forward = path.join(Config.HYSPLIT_folder, sourcename, "HYSPLIT_GFSQ_F_2019091423_2019091504.txt")
    if not logger is None: logger.info("Processing file {}".format(t_backward))
    ah.addTrajectories(Algorithm_HYSPLIT_from_trajectory(source, model, chart_datetime, t_backward).trajectories)
    ah.addTrajectories(Algorithm_HYSPLIT_from_trajectory(source, model, chart_datetime, t_forward).trajectories)

    if 0 == len(ah.trajectories):
        if not logger is None: logger.error("No HYSPLIT trajectory files found")
        print("No HYSPLIT trajectory files found")

    if not fileout is None:
        ah.chart_time_series(Geometry.CONTAINS, orbit, processor, fileout)
    else:
        ah.chart_time_series(Geometry.CONTAINS, orbit, processor, None)

def _chart_IME(sourcename: str, orbit: str, processor: str, background: float, transectid: int, imagefile: str) -> None:
    filename = Algorithm_CSF.get_picklename(Config(), sourcename, orbit, processor)
    
    if not path.exists(filename):
        __utility_file_does_not_exist(filename)
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    ime = Algorithm_IME_factory(sad, background, transectid)
    ime.chart(imagefile)

def _chart_sadavarte_figure2(type:str, fileout:str) -> int:
    ''' 
    Reproduce Figure 2 of Sadaverte 2021 using data from local analysis 
    For all filtered orbits, 
    1. open a pickle file, 
    2. check is status is success
    3. check if configuration is correct, 
    4. add to dates_string, add to counts
    
    Parameters
    ----------
    type: str
        - "figure2_all" - no filter, 
        - "figure2_min_max" - additionally remove fixed large orbits
    '''
    with open(TROPOMI_Filter.pickle, "rb") as f:
        filter: TROPOMI_Filter = load(f)
        f.close()

    config = Config()
    orbits = filter.result.copy()
    ncfiles = filter.find_highest_configured_procesor(config.TROPOMI_Processors)

    dates = []
    emissionrates = []
    config = Config()
    for i in orbits:
        strorbit = str(i).zfill(5)
        processor = TROPOMI.get_processor_version(ncfiles[i])
        pickle_run_filename_full = Algorithm_CSF.get_picklename(config, "HailCreek", strorbit, processor)
        
        if not path.exists(pickle_run_filename_full):
            if not logger is None: logger.error("Missing run data for orbit: {}".format(strorbit))
            continue

        with open(pickle_run_filename_full, "rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()

        if not sad.status == Status.SUCCESS.value:
            if not logger is None: logger.error("Algorithm run data for orbit: {} has not completed succesfully. Status {} ".format(i, sad.status))
            continue

        if type == "figure2_min_max":
            if sad.transects_valid_positive_pixels_count < config.Algorithm_CSF_transect_positive_minimum_count: continue
            if config.Algorithm_CSF_transect_valid_maximumu_count < sad.transects_valid_pixels_count: continue
            # if 0 < sad.count_source_in_downindbox(): continue
            # if 0 < sad.count_source_in_upwindbox(): continue
        
        if not logger is None: logger.info("Figure 2: Orbit included: {} Date: {}".format(sad.tropomi.orbit, sad.tropomi_source_date))
        print("Figure 2: Orbit included: {} Date: {}".format(sad.tropomi.orbit, sad.tropomi_source_date))

        dates.append(sad.tropomi_source_date)
        emissionrates.append(sad.q_valid_hour) 

    count = len(emissionrates)
    avg = average(emissionrates)
    stddev = std(emissionrates)
    sem2 = 2 * stddev / sqrt(count) # Symmetrical confidence interval for standar error of mean (sem) at 95% confidence level
    xticks = [
        datetime.strptime("2018-06-01", "%Y-%m-%d"),
        datetime.strptime("2018-09-01", "%Y-%m-%d"),
        datetime.strptime("2018-12-01", "%Y-%m-%d"),
        datetime.strptime("2019-03-01", "%Y-%m-%d"),
        datetime.strptime("2019-06-01", "%Y-%m-%d"),
        datetime.strptime("2019-09-01", "%Y-%m-%d"),
        datetime.strptime("2019-12-01", "%Y-%m-%d")
    ]

    xlabels = [
        "Jun\n2018",
        "Sep\n2018",
        "Dec\n2018",
        "Mar\n2019",
        "Jun\n2019",
        "Sep\n2019",
        "Dec\n2019"
    ]

    _, ax = plt.subplots(figsize=(10,5.5), layout=None)
    ax.scatter(dates, emissionrates)
    ax.set_xticks(xticks, labels=xlabels)
    ax.set_yticks(ticks=[-20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110], labels=["-20", "", "0", "", "20", "", "40", "", "60", "", "80", "", "100", ""])    
    ax.set_xlabel("Date", loc="center")
    ax.set_ylabel("Methane flux (t hr-1)", loc="center")
    
    txtstr = "\n".join((
        "Count: {}".format(count), 
        r"$\mu: {:.3f} (t/h)$".format(avg),
        r"$2 * \sigma: {:.3f} (t/h)$".format(sem2))
    )
    ax.text(0.01,0.95, txtstr, transform=ax.transAxes, verticalalignment="top", bbox=dict(boxstyle="square", facecolor="wheat", alpha=0.5))

    if not fileout is None:
        plt.savefig(fileout )
        print("Chart generated: {}".format(fileout))
    else:
        plt.show()
    plt.close("all")
    return

def _chart_SRTM(t_lat: float, l_lon: float, b_lat:float, r_lon:float, fileout:str) -> None:
    '''
    Chart SRTM data in the box

    Parameters
    ----------
    t_lat: float
        Top latitude
    l_lon: float
        Left longitude
    b_lat:float
        Bottom latitude
    r_lon:float
        Right longitude
    fileout: path
        If specified location of output

    '''
    if not path.exists(path_9s_geotiff):
        __utility_file_does_not_exist(path_9s_geotiff)
        return 
    
    poly = box(l_lon, b_lat, r_lon, t_lat)
    srtm = SRTM_from_geotiff(path_9s_geotiff, t_lat, l_lon, b_lat, r_lon)
    s_h = create_source("HailCreek")
    s_mn = create_source("MoranbahNorth")
    aws = create_aws("Moranbah")

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())

    srtm.chart(fig, ax)
    ax.coastlines()
    # Chart sources that have activity data between April 2018 and December 2020
    for s in Source.Sources:
        source = create_source(s)
        if (not 2018 in source.activity) and (not 2019 in source.activity) and (not 2020 in source.activity): continue
        if poly.contains(source.xy): source.plot_source(ax)
    # if poly.contains(s_h.xy): s_h.plot_source(ax)
    # if poly.contains(s_mn.xy): s_mn.plot_source(ax)
    if poly.contains(aws.xy): aws.plot_aws(ax)

    plt.title("SRTM 9s Topography (m) Hail Creek")

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.left_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'gray'}
    gl.ylabel_style = {'size': 10, 'color': 'grey'} 

    if not fileout is None:        
        plt.savefig(fileout)
        print("Chart generated: {}".format(fileout))
    else:
        plt.show()

    plt.close("all")
    return 

def _chart_synthetic_geometry(fileout: str) -> None:
    '''
    Draw chart illustrating geometries
    '''
    fig, ax = plt.subplots(figsize=(10,5))
    plt.title("Geometrical relationships: CONTAINS, CENTER, INTERSECTS")
    ax.set_xlim(xmin=0, xmax = 10)
    ax.set_ylim(ymin = 0, ymax = 5)

    c_x = 5.0
    c_y = 2.5

    v_1_x = [2.0, 2.0]
    v_1_y = [1.5, 3.5]

    v_2_x = [2.0 + 7.0 / 5.0 , 2.0 + 7.0 / 5.0 ]
    v_2_y = [1.5, 3.5]

    m_1_x = (v_1_x[0] + v_2_x[0]) / 2.0
    m_1_y = (v_1_y[0] + v_1_y[1]) / 2.0

    v_3_x = [2.0 + 14.0 / 5.0, 2.0 + 14.0 / 5.0]
    v_3_y = [1.5, 3.5]

    m_2_x = (v_2_x[0] + v_3_x[0]) / 2.0
    m_2_y = (v_2_y[0] + v_2_y[1]) / 2.0

    v_4_x = [2.0 + 21.0 / 5.0, 2.0 + 21.0 / 5.0]
    v_4_y = [1.5, 3.5]

    m_3_x = (v_3_x[0] + v_4_x[0]) / 2.0
    m_3_y = (v_3_y[0] + v_3_y[1]) / 2.0

    h_1_x = [1.5, 7.0]
    h_1_y = [2.0, 2.0]

    h_2_x = [1.5, 7.0]
    h_2_y = [3.0, 3.0]

    a_1_x = [v_1_x[0] - 0.5, v_1_x[0] - 0.5]
    a_1_y = [h_1_y[0] + 0.2, h_2_y[0] - 0.2]

    theta = 17 * pi / 180

    # Rotate points of the first vertical line
    a_1_x[0], a_1_y[0] = __utility_rotate(a_1_x[0], a_1_y[0], c_x, c_y, theta)
    a_1_x[1], a_1_y[1] = __utility_rotate(a_1_x[1], a_1_y[1], c_x, c_y, theta)

    v_1_x[0], v_1_y[0] = __utility_rotate(v_1_x[0], v_1_y[0], c_x, c_y, theta)
    v_1_x[1], v_1_y[1] = __utility_rotate(v_1_x[1], v_1_y[1], c_x, c_y, theta)

    v_2_x[0], v_2_y[0] = __utility_rotate(v_2_x[0], v_2_y[0], c_x, c_y, theta)
    v_2_x[1], v_2_y[1] = __utility_rotate(v_2_x[1], v_2_y[1], c_x, c_y, theta)

    v_3_x[0], v_3_y[0] = __utility_rotate(v_3_x[0], v_3_y[0], c_x, c_y, theta)
    v_3_x[1], v_3_y[1] = __utility_rotate(v_3_x[1], v_3_y[1], c_x, c_y, theta)

    v_4_x[0], v_4_y[0] = __utility_rotate(v_4_x[0], v_4_y[0], c_x, c_y, theta)
    v_4_x[1], v_4_y[1] = __utility_rotate(v_4_x[1], v_4_y[1], c_x, c_y, theta)

    h_1_x[0], h_1_y[0] = __utility_rotate(h_1_x[0], h_1_y[0], c_x, c_y, theta)
    h_1_x[1], h_1_y[1] = __utility_rotate(h_1_x[1], h_1_y[1], c_x, c_y, theta)

    h_2_x[0], h_2_y[0] = __utility_rotate(h_2_x[0], h_2_y[0], c_x, c_y, theta)
    h_2_x[1], h_2_y[1] = __utility_rotate(h_2_x[1], h_2_y[1], c_x, c_y, theta)

    m_1_x, m_1_y = __utility_rotate(m_1_x, m_1_y, c_x, c_y, theta)
    m_2_x, m_2_y = __utility_rotate(m_2_x, m_2_y, c_x, c_y, theta)
    m_3_x, m_3_y = __utility_rotate(m_3_x, m_3_y, c_x, c_y, theta)

    rect = patches.Rectangle((1.0, 1.1), 6, 1.5, linewidth=1, edgecolor='brown', facecolor='none')
    h_1 = lines.Line2D(h_1_x, h_1_y, linewidth=1, color="blue")
    h_2 = lines.Line2D(h_2_x, h_2_y, linewidth=2, color="blue")
    v_1 = lines.Line2D(v_1_x, v_1_y, linewidth=2, color="blue")
    v_2 = lines.Line2D(v_2_x, v_2_y, linewidth=2, color="blue")
    v_3 = lines.Line2D(v_3_x, v_3_y, linewidth=2, color="blue")
    v_4 = lines.Line2D(v_4_x, v_4_y, linewidth=2, color="blue")

    ax.set_axis_off()
    ax.add_patch(rect)
    ax.add_line(h_1)
    ax.add_line(h_2)
    ax.add_line(v_1)
    ax.add_line(v_2)
    ax.add_line(v_3)
    ax.add_line(v_4)
    ax.arrow(a_1_x[0], a_1_y[0], dx=a_1_x[1] - a_1_x[0], dy = a_1_y[1] - a_1_y[0], head_width=0.1, head_length=0.1, length_includes_head=True)
    ax.scatter([m_1_x, m_2_x, m_3_x], [m_1_y, m_2_y, m_3_y])
    ax.text(m_1_x, m_1_y + 0.1,"(i,j)")
    ax.text(m_2_x, m_2_y + 0.1,"(i,j+1)")
    ax.text(m_3_x, m_3_y + 0.1,"(i,j+2)")
    ax.text(a_1_x[1] + 0.15, a_1_y[1] - 0.25, "Flight\ndirection", rotation=17-90, rotation_mode="anchor", size=9, horizontalalignment="center", verticalalignment="bottom")

    if not fileout is None:
        plt.savefig(fileout)
        print("Chart generated: {}".format(fileout))
    else:
        plt.show()
    
    plt.close("all")
    return

def _chart_synthetic_masking(fileout: str) -> None:
    '''
    Draw chart illustrating masking
    '''
    lcx = 0.1
    lcy = 0.25
    rw = 0.3
    rh = 0.25
    txt_offset_x = rw * 0.5
    txt_offset_y = rh * 0.5
    
    b_x = [lcx, lcx + rw, lcx + rw + rw]
    b_y = [lcy, lcy + rh]

    # Top row left chart
    rect_l_t_0 = patches.Rectangle((b_x[0], b_y[1]), rw, rh, linewidth=1, edgecolor='white', facecolor='blue')
    rect_l_t_1 = patches.Rectangle((b_x[1], b_y[1]), rw, rh, linewidth=1, edgecolor='white', facecolor='blue')
    rect_l_t_2 = patches.Rectangle((b_x[2], b_y[1]), rw, rh, linewidth=1, edgecolor='white', facecolor='blue')
    
    # Bottom row left chart
    rect_l_b_0 = patches.Rectangle((b_x[0], b_y[0]), rw, rh, linewidth=1, edgecolor='white', facecolor='brown')
    rect_l_b_1 = patches.Rectangle((b_x[1], b_y[0]), rw, rh, linewidth=1, edgecolor='white', facecolor='brown')
    rect_l_b_2 = patches.Rectangle((b_x[2], b_y[0]), rw, rh, linewidth=1, edgecolor='white', facecolor='brown')

    # Top row right chart
    rect_r_t_0 = patches.Rectangle((b_x[0], b_y[1]), rw, rh, linewidth=1, edgecolor='white', facecolor='green')
    rect_r_t_1 = patches.Rectangle((b_x[1], b_y[1]), rw, rh, linewidth=1, edgecolor='white', facecolor='blue')
    rect_r_t_2 = patches.Rectangle((b_x[2], b_y[1]), rw, rh, linewidth=1, edgecolor='white', facecolor='blue')
    
    # Bottom row right chart
    rect_r_b_0 = patches.Rectangle((b_x[0], b_y[0]), rw, rh, linewidth=1, edgecolor='white', facecolor='brown')
    rect_r_b_1 = patches.Rectangle((b_x[1], b_y[0]), rw, rh, linewidth=1, edgecolor='white', facecolor='brown')
    rect_r_b_2 = patches.Rectangle((b_x[2], b_y[0]), rw, rh, linewidth=1, edgecolor='white', facecolor='brown')

    fig, axs = plt.subplots(figsize=(10,3), nrows=1, ncols=2)
    fig.suptitle("Masking: NONE, NEGATIVE")

    ax0 = axs[0]
    ax0.set_axis_off()
    ax0.add_patch(rect_l_b_0)
    ax0.text(b_x[0] + txt_offset_x, b_y[0] + txt_offset_y, "10 ppb", color="white", horizontalalignment="center", verticalalignment="center" )
    ax0.add_patch(rect_l_b_1)
    ax0.text(b_x[1] + txt_offset_x, b_y[0] + txt_offset_y, "10 ppb", color="white", horizontalalignment="center", verticalalignment="center" )
    ax0.add_patch(rect_l_b_2)
    ax0.text(b_x[2] + txt_offset_x, b_y[0] + txt_offset_y, "10 ppb", color="white", horizontalalignment="center", verticalalignment="center" )

    ax0.add_patch(rect_l_t_0)
    ax0.text(b_x[0] + txt_offset_x, b_y[1] + txt_offset_y, "-3 ppb", color="white", horizontalalignment="center", verticalalignment="center" )
    ax0.add_patch(rect_l_t_1)
    ax0.text(b_x[1] + txt_offset_x, b_y[1] + txt_offset_y, "3 ppb", color="white", horizontalalignment="center", verticalalignment="center" )
    ax0.add_patch(rect_l_t_2)
    ax0.text(b_x[2] + txt_offset_x, b_y[1] + txt_offset_y, "4 ppb", color="white", horizontalalignment="center", verticalalignment="center" )

    ax1 = axs[1]
    ax1.set_axis_off()
    ax1.add_patch(rect_r_b_0)
    ax1.text(b_x[0] + txt_offset_x, b_y[0] + txt_offset_y, "10 ppb", color="white", horizontalalignment="center", verticalalignment="center" )
    ax1.add_patch(rect_r_b_1)
    ax1.text(b_x[1] + txt_offset_x, b_y[0] + txt_offset_y, "10 ppb", color="white", horizontalalignment="center", verticalalignment="center" )
    ax1.add_patch(rect_r_b_2)
    ax1.text(b_x[2] + txt_offset_x, b_y[0] + txt_offset_y, "10 ppb", color="white", horizontalalignment="center", verticalalignment="center" )

    ax1.add_patch(rect_r_t_0)
    ax1.text(b_x[0] + txt_offset_x, b_y[1] + txt_offset_y, "-3 ppb", color="white", horizontalalignment="center", verticalalignment="center" )
    ax1.add_patch(rect_r_t_1)
    ax1.text(b_x[1] + txt_offset_x, b_y[1] + txt_offset_y, "3 ppb", color="white", horizontalalignment="center", verticalalignment="center" )
    ax1.add_patch(rect_r_t_2)
    ax1.text(b_x[2] + txt_offset_x, b_y[1] + txt_offset_y, "4 ppb", color="white", horizontalalignment="center", verticalalignment="center" )

    if not fileout is None:
        plt.savefig(fileout)
        print("Chart generated: {}".format(fileout))
    else:
        plt.show()
    
    plt.close("all")
    return

def _chart_synthetic_transect(fileout: str) -> None:
    '''
    Draw chart of synthetic data to demonstrate impact of selection of background value
    '''
    # Concentration function
    # c = (-1. * (t - 1.0) * (t -5.0) * (t - 7.0) * (t - 10.0) / 40.) + 1830.

    # Coefficients
    scale = 40.
    c4 = -1.
    c3 = 1. + 5. + 7. + 10.
    c2 = -( 1. * 5. + 1. * 7. + 1. * 10. + 5. * 7. + 5. * 10. + 7. * 10.)
    c1 = ( 1. * 5. * 7. + 1. * 5. * 10. +  1. * 7. * 10. + 5. * 7. * 10.)
    c0 = 1830.

    # Concentration function numerical values
    t = linspace(1, 10, 100)
    
    c = (c4 * pow(t, 4) + c3 * pow(t,3) + c2 * pow(t, 2) + c1 * t) / scale + c0

    # print("{}x^4 + {}x^3+ {}x^2 + {}x +{}".format(c4, c3, c2, c1, c0))

    # Integrated emissions
    barx = linspace(1, 10, 10)
    integral = []
    for i in range(10):
        x_s = float(i + 0.5)
        x_e = float(i + 1.5)
        i_s = (c4 * pow(x_s, 5) / 5. + c3 * pow(x_s, 4) / 4. + c2 * pow(x_s, 3) / 3. + c1 * pow(x_s, 2) / 2.) / scale + c0 * x_s
        i_e = (c4 * pow(x_e, 5) / 5. + c3 * pow(x_e, 4) / 4. + c2 * pow(x_e, 3) / 3. + c1 * pow(x_e, 2) / 2.) / scale + c0 * x_e
        integral.append(i_e - i_s)

    ymin = min(integral)
    ymax = max(integral)
    # print("min: {}, max: {}".format(ymin, ymax))

    # Possible backgrounds
    b0 = 1841.0
    b1 = 1839.5
    b2 = 1838.0
    e0 = 0.
    e1 = 0.
    e2 = 0.
    for i in range(10):
        if integral[i] > b0 : e0 += integral[i] - b0
        if integral[i] > b1 : e1 += integral[i] - b1
        if integral[i] > b2 : e2 += integral[i] - b2

    style = {'facecolor': 'none', 'edgecolor': 'C0', 'linewidth': 3}

    fig = plt.figure(figsize=(10,10))
    gs = GridSpec(4, 1, height_ratios=[5, 1, 1, 1], wspace=0.05)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    ax4 = fig.add_subplot(gs[3])

    fig.suptitle("Synthetic concentration, satellite pixel values\n" +
    "Pixels in plume transect for background: {:4.0f} (ppb) => Enhancement: {:.3f} (ppb)\n".format(b0, e0) +
    "Pixels in plume transect for background: {:4.0f} (ppb) => Enhancement: {:.3f} (ppb)\n".format(b1, e1) +
    "Pixels in plume transect for background: {:4.0f} (ppb) => Enhancement: {:.3f} (ppb)".format(b2, e2))
    ax1.set_xlim(xmin=1, xmax = 10)
    ax1.set_ylim(ymin = ymin - 1, ymax = ymax + 1)
    ax1.plot(t, c, linewidth = 2, color="green")
    ax1.bar(barx, integral, width=1.0, **style )
    ax1.hlines(b0, xmin=1, xmax = 10, linewidth = 2, color="pink")
    ax1.hlines(b1, xmin=1, xmax = 10, linewidth = 2, color="red")
    ax1.hlines(b2, xmin=1, xmax = 10, linewidth = 2, color="brown")
    ax1.legend(["Concentration", "Background 1841 (ppb)", "Background 1839.5 (ppb)", "Background 1838 (ppb)", "Satellite pixel"])

    # Shape of pixels intersection
    v2 =  [1,1,1,1,1,1,1,1,1,1]
    c2 =  ['white', 'pink', 'pink', 'white','white', 'white', 'white', 'white','white', 'white']
    style2 = {'edgecolor': 'C0', 'linewidth': 3}
    ax2.set_xlim(xmin=1, xmax = 10)
    ax2.set_ylim(ymin=0, ymax=1)
    ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax2.bar(barx, v2, color=c2, width=1.0, **style2)

    v3 =  [1,1,1,1,1,1,1,1,1,1]
    c3 =  ['white', 'red', 'red', 'red','white', 'white', 'white', 'red','red', 'white']
    style3 = {'edgecolor': 'C0', 'linewidth': 3}
    ax3.set_xlim(xmin=1, xmax = 10)
    ax3.set_ylim(ymin=0, ymax=1)
    ax3.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax3.bar(barx, v3, color=c3, width=1.0, **style3)

    v4 =  [1,1,1,1,1,1,1,1,1,1]
    c4 =  ['brown', 'brown', 'brown', 'brown','brown', 'brown', 'brown', 'brown','brown', 'brown']
    style4 = {'edgecolor': 'C0', 'linewidth': 3}
    ax4.set_xlim(xmin=1, xmax = 10)
    ax4.set_ylim(ymin=0, ymax=1)
    ax4.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax4.bar(barx, v4, color=c4, width=1.0, **style4)

    if not fileout is None:
        plt.savefig(fileout)
        print("Chart generated: {}".format(fileout))
    else:
        plt.show()
    
    plt.close("all")

def _chart_TM(sourcename: str, orbit: str, processor: str, model: str, direction: str, background: float, imagefile: str) -> None:
    source = create_source(sourcename)
    filename = Algorithm_CSF.get_picklename(Config(), sourcename, orbit, processor)
    
    if not path.exists(filename):
        __utility_file_does_not_exist(filename)
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    sourcepath = path.join(Config.HYSPLIT_folder, source.case_name)
    filenames = listdir(sourcepath)
    pattern = "HYSPLIT_{}_{}_".format(model, direction)
    date = datetime(sad.tropomi_source_date.year, sad.tropomi_source_date.month, sad.tropomi_source_date.day, sad.tropomi_source_date.hour, tzinfo=timezone.utc)
    ah = Algorithm_HYSPLIT(source, model, date)
    for f in filenames:
        if not f.startswith(pattern): continue
        if not logger is None: logger.info("Processing file {}".format(f))
        ah.addTrajectories(Algorithm_HYSPLIT_from_trajectory(source, model, date, path.join(sourcepath, f)).trajectories)

    if 0 == len(ah.trajectories):
        if not logger is None: logger.error("No HYSPLIT trajectory files found")
        print("No HYSPLIT trajectory files found")
        return
    
    tm = Algorithm_TM_factory(sad, ah, background)
    tm.chart(imagefile)
    return

def _list_activity_HailCreek_CY(fileout: str) -> None:
    '''
    List activity data coming from Qld Statistics. 2025-calendar-year-coal-production-statistics.xls
    '''
    __utility_table_output_csv(fileout, csv_activity_cy)
    return

def _list_activity_HailCreek_FY(fileout: str) -> None:
    '''
    List activity data coming from Qld Statistics. coal-production-data-fy2010-fy2025.xlsx
    '''
    __utility_table_output_csv(fileout, csv_activity_fy)
    return

def _list_BOM_Moranbabah(date_start: datetime, date_end: datetime, fileout: str) -> None:
    '''
    Generate a CSV file of Moranbah AWS data for specified period
    Date, Air T C(deg), Dew Point C(deg), Wind Speed (km/h), Wind direction in degrees true
    '''
    if not path.exists(path_Moranbah_AWS_CSV):
        __utility_file_does_not_exist(path_Moranbah_AWS_CSV)
        return

    aws = BOM_AWS_from_csv("Moranbah AWS", path_Moranbah_AWS_CSV)
    csv = aws.list_BOM_Moranbabah(date_start, date_end)
    __utility_table_output_csv(fileout, csv)

    return

def _list_ERA5_BLH(source_casename: str, orbit:str, processor: str):
    '''
    List boundary layer height at the ERA5 grid point closest to HailCreek mine at 04Z
    '''
    tropomi = TROPOMI_for_orbit(orbit, processor)
    source = create_source(source_casename)
    era5_sl = ERA5_SingleLevel.from_netCDF(ERA5_SingleLevel.get_filepath("CSF", source.case_name, orbit))
    [_, _, dt, grid_timestamp] = tropomi.get_pixel_for_source(source.xy)
    grid_time = ERA5.calculate_ERA_GridTime(dt)
    era5_sl.convert_timestamp_to_index(grid_time)
    era5_sl.convert_lat_to_index(source.gridpoint.y)
    era5_sl.convert_lon_to_index(source.gridpoint.x)

    blh = era5_sl.get("blh",era5_sl.time_index, era5_sl.lat_index, era5_sl.lon_index)
    
    tstring = grid_timestamp.strftime("%Y %m %d %H:%M UTC")
    print("Time: {}, Latitude: {}, Longitude: {}".format(tstring, source.gridpoint.y, source.gridpoint.x))
    print("BLH (m): {}".format(blh))

def _list_ERA5_SP(source_casename: str, orbit:str, processor: str):
    '''
    List surface pressure at the ERA5 grid point closest to HailCreek mine at 04Z
    '''
    tropomi = TROPOMI_for_orbit(orbit, processor)
    source = create_source(source_casename)
    era5_sl = ERA5_SingleLevel.from_netCDF(ERA5_SingleLevel.get_filepath("CSF", source.case_name, orbit))
    [_, _, dt, grid_timestamp] = tropomi.get_pixel_for_source(source.xy)
    grid_time = ERA5.calculate_ERA_GridTime(dt)
    era5_sl.convert_timestamp_to_index(grid_time)
    era5_sl.convert_lat_to_index(source.gridpoint.y)
    era5_sl.convert_lon_to_index(source.gridpoint.x)
    sp = era5_sl.get("sp",era5_sl.time_index, era5_sl.lat_index, era5_sl.lon_index)

    tstring = grid_timestamp.strftime("%Y %m %d %H:%M UTC")
    print("Time: {}, Latitude: {}, Longitude: {}".format(tstring, source.gridpoint.y, source.gridpoint.x))
    print("Surface Pressure (Pa): {}".format(sp))

def _list_ERA5_Wind(source_casename: str, orbit:str, processor: str, fileout: str):
    '''
    List wind user for  at the ERA5 grid point closest to HailCreek mine at 04Z
    '''
    tropomi = TROPOMI_for_orbit(orbit, processor)
    source = create_source(source_casename)
    era5_p = ERA5_PressureLevels.from_netCDF(ERA5_PressureLevels.get_filepath("CSF",source.case_name, orbit))
    [_,_,dt, _] = tropomi.get_pixel_for_source(source.xy)
    grid_time = ERA5.calculate_ERA_GridTime(dt)
    grid_point: Point = ERA5.calculate_ERA5_NearestGridPoint(source.xy)

    # Find vertical wind profile om 20190915 04Z at point nearest to mine
    era5_p.convert_timestamp_to_index(grid_time)
    era5_p.convert_lat_to_index(grid_point.y)
    era5_p.convert_lon_to_index(grid_point.x)

    # Pressure, geopotential, u and v as 1-d arrays over pressure levels
    p_wind = era5_p.pressure
    z_wind = era5_p.variables["z"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]
    u = era5_p.variables["u"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]
    v = era5_p.variables["v"][era5_p.time_index, :, era5_p.lat_index, era5_p.lon_index]

    pl = len(p_wind)
    csv = ["Pressure (hpa), Height (m), U (m/s), V (m/s)"]
    for p in range(pl): csv.append("{}, {}, {}, {}".format(p_wind[p], z_wind[p]/GRAVITY, u[p], v[p]))

    if fileout is None:
        tstring = datetime.fromtimestamp(grid_time).strftime("%Y %m %d %H:%M UTC")
        print("Time: {}, Latitude: {}, Longitude: {}".format(tstring, grid_point.y, grid_point.x))

    __utility_table_output_csv(fileout, csv)

def _list_csf_tm_comparison(sourcename: str, orbit: str, processor: str, minbackground: float, maxbackground: float, fileout: str) -> None:
    '''
    Run comparison of cross sectional flux and integrated mass enhancement method for orbit
    '''
    backgrounds = range(int(minbackground), int(maxbackground))    
    model = "GFSQ"
    direction = "F"

    # Get Algorithm_CSF data
    filename = Algorithm_CSF.get_picklename(Config(), sourcename, orbit, processor)
    if not path.exists(filename):
        __utility_file_does_not_exist(filename)
        return

    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    # Create HYSPLIT
    source = create_source("HailCreek")
    sourcepath = path.join(Config.HYSPLIT_folder, source.case_name)
    filenames = listdir(sourcepath)
    pattern = "HYSPLIT_{}_{}_".format(model, direction)
    date = datetime(sad.tropomi_source_date.year, sad.tropomi_source_date.month, sad.tropomi_source_date.day, sad.tropomi_source_date.hour, tzinfo=timezone.utc)
    ah = Algorithm_HYSPLIT(source, model, date)

    for f in filenames:
        if not f.startswith(pattern): continue
        if not logger is None: logger.info("Processing file {}".format(f))
        ah.addTrajectories(Algorithm_HYSPLIT_from_trajectory(source, model, date, path.join(sourcepath, f)).trajectories)

    if 0 == len(ah.trajectories):
        m = "Error: No HYSPLIT trajectory files found"
        print(m)
        if not logger is None: logger.error(m)
        return
    
    csv=["Background (ppb), CSF (t/h), TM (t/h)"]

    for b in backgrounds:
        sad.status = Status.PROCESSING.value
        sad.background = b
        sad.calculate_transects()
        sad.calculate_sourcerate()
        tm = Algorithm_TM_factory(sad, ah, None) # Use background set is sad
        if sad.status == Status.SUCCESS.value:
            csv.append("{},{},{}".format(b, sad.q_valid_hour, tm.q_valid_hour))
        else:
            csv.append("{},ERROR,{}".format(b, tm.q_valid_hour))
    
    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_Figure2(type: str, fileout:str) -> None:
    ''' 
    Reproduce Figure 2 of Sadaverte 2021 using data from local analysis 
    For all filtered orbits, 
    1. open a pickle file, 
    2. check is status is success
    3. check if configuration is correct, 
    4. add to dates_string, add to counts

    Parameters
    ----------
    type: str
        - "figure2_all" - no filter, 
        - "figure2_min_max" - additionally remove fixed large orbits
    '''
    if not path.exists(TROPOMI_Filter.pickle): 
        __utility_file_does_not_exist(TROPOMI_Filter.pickle)
        return
    
    with open(TROPOMI_Filter.pickle, "rb") as f:
        filter: TROPOMI_Filter = load(f)
        f.close()

    config = Config()
    orbits = filter.result.copy()
    ncfiles = filter.find_highest_configured_procesor(config.TROPOMI_Processors)
    source = create_source("HailCreek")

    failed = []
    csv = ["orbit, date, q (t/hr), transect valid pixel count, transect positive pixel count"]
    values = []
    config = Config()
    count_passed = 0
    count_failed = 0
    for i in orbits:
        strorbit = str(i).zfill(5)
        processor = TROPOMI.get_processor_version(ncfiles[i])
        pickle_run_filename_full = Algorithm_CSF.get_picklename(config, "HailCreek", strorbit, processor)
        
        if not path.exists(pickle_run_filename_full):
            failed.append("{} Missing run data".format(strorbit))
            continue

        with open(pickle_run_filename_full, "rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()

        if not sad.status == Status.SUCCESS.value:
            failed.append("{:05d} Algorithm failed with status {} ".format(i, sad.status))
            count_failed += 1
            continue

        if  type == "figure2_min_max":
            if sad.transects_valid_positive_pixels_count < config.Algorithm_CSF_transect_positive_minimum_count:
                failed.append("{:05d} Transect positive pixel count {} < {}. Orbit skipped.".format(i, sad.transects_valid_positive_pixels_count, config.Algorithm_CSF_transect_positive_minimum_count))
                count_failed += 1
                continue

            if config.Algorithm_CSF_transect_valid_maximumu_count < sad.transects_valid_pixels_count:
                failed.append("{:05d} {} < Transact valid pixel count {}. Orbit skipped.".format(i, config.Algorithm_CSF_transect_valid_maximumu_count, sad.transects_valid_pixels_count))
                count_failed += 1
                continue

            # if 0 < sad.count_source_in_downindbox():
            #     failed.append("{:05d} Sources present in the downwind box. Orbit skipped.".format(i))
            #     count_failed += 1
            #     continue

            # if 0 < sad.count_source_in_upwindbox():
            #     failed.append("{:05d} Sources present in the upwind box. Orbit skipped.".format(i))
            #     count_failed += 1
            #     continue

        v = sad.q_valid_hour
        csv.append("{:05d}, {}-{:02d}-{:02d}, {}, {}, {}".format(i, sad.tropomi_source_date.year, sad.tropomi_source_date.month, sad.tropomi_source_date.day, v, sad.transects_valid_pixels_count, sad.transects_valid_positive_pixels_count))
        values.append(v)
        count_passed += 1

    q_h = average(values)
    std_h2 = 2 * std(values)/ sqrt(count_passed) # standard error of the mean this corresponds to 95% confidence interval of normal distribution
    q_y = q_h * 24 * 365 * 1E-3
    sem_y2 = std_h2 * 24 * 365 * 1E-3
    # Note that 2019 is from July 2018 - June 2019, 
    ief = (2 * q_y) / (source.activity[2019] + source.activity[2020])

    if fileout is None:
        messages = []
        if type == "figure2_all": messages.append("All: ")
        if type == "figure2_min_max": messages.append("Min and Max: ")

        messages.append("Count total: {} Count Failed: {} Count Passed: {}".format(count_passed + count_failed, count_failed, count_passed))
        messages.append("Orbit average q: {:.3f} (t/hour) 2 * sem: {:.3f} (t/h)".format(q_h, std_h2))
        messages.append("Annual average annual q: {:.3f} (Gg/year) 2 * sem: {:.3f}".format(q_y, sem_y2))
        messages.append("Activity 2019: {:.3f} (Mt) Activity 2020: {:.3f} (Mt) ief: {:.3f} (g/kg)".format(source.activity[2019], source.activity[2020], ief))
        
        for m in failed: print(m)
        print()

    __utility_table_output_csv(fileout, csv)

    if fileout is None:
        for m in messages: print(m)
        print()

    return

def _list_HailCreek(fileout: str) -> None:
    '''
    List basic data about Hail Creek
    '''
    __utility_table_output_csv(fileout, csv_hailcreek)
    return

def _list_HailCreek_emissions(fileout: str) -> None:
    '''
    List emissions Hail Creek
    '''
    __utility_table_output_csv(fileout, csv_emissions)
    return

def _list_tropomi_filteredorbits(type: str, fileout: str) -> None:
    '''
    List content of filters into csv files
    '''

    if not path.exists(TROPOMI_Filter.pickle):
        m = "Error: Filter file: {} is missing. Please run TROPOMI_Filter program first.".format(TROPOMI_Filter.pickle)
        print(m)
        if not logger is None: logger.error(m)
        return

    with open(TROPOMI_Filter.pickle, "rb") as f:
        filter: TROPOMI_Filter = load(f)
        f.close()

    csv: list[str] = ["Orbit, File"]
    lt = len(filter.success_contains)
    for i in range(lt):
        orbit = filter.success_contains[i]
        file = filter.success_contains_s3_key[i].replace(".xml", ".nc")
        if type == "contains": 
            csv.append("{},{}".format(orbit, file))       
        elif type == "qa" and orbit in filter.success_qa:
            csv.append("{},{}".format(orbit, file))
        elif type == "r" and orbit in filter.success_r:
            csv.append("{},{}".format(orbit, file))  

    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_comparison_background(source:str, orbit:str, processor: str, fileout:str) -> None:
    ''' 
    Generate a csv file and print out comparison of background for different configurations
    '''
    models = []
    config111 = Config()
    config111.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    models.append(Algorithm_CSF.get_picklename(config111, source, orbit, processor))

    config211 = Config()
    config211.Algorithm_CSF_background_geometry = Geometry.CENTER
    models.append(Algorithm_CSF.get_picklename(config211, source, orbit, processor))

    config311 = Config()
    config311.Algorithm_CSF_background_geometry = Geometry.INTERSECTS
    models.append(Algorithm_CSF.get_picklename(config311, source, orbit, processor))
    
    csv = []
    csv.append("Background geometry, Pixel count, Average (ppb)")
    for m in models:
        if not path.exists(m):
            message = "Model {} is missing. Please run models first.".format(m)
            print(message)
            if not logger is None: logger.error(message)
            return

        with open(m, "rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()

        csv.append("{}, {}, {}".format(sad.config.Algorithm_CSF_background_geometry.name, sad.upwind_box_background_count, sad.upwind_box_background_average))

    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_comparison_emissionrate(source:str, orbit: str, processor: str, fileout:str) -> None:
    ''' 
    Generate a csv file and print out comparison of emission reates for different configurations

    Parameters
    ----------
    source: str
        Source name
    orbit: str
        Orbit for analysis
    '''
    models = []

    config111 = Config() 
    config111.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config111.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config111.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config111, source, orbit, processor))
    
    config112 = Config()
    config112.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config112.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config112.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config112, source, orbit, processor))

    config121 = Config() 
    config121.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config121.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config121.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config121, source, orbit, processor))

    config122 = Config()
    config122.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config122.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config122.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config122, source, orbit, processor))

    config131 = Config() 
    config131.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config131.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config131.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config131, source, orbit, processor))
    
    config132 = Config()
    config132.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config132.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config132.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config132, source, orbit, processor))

    config211 = Config() 
    config211.Algorithm_CSF_background_geometry = Geometry.CENTER
    config211.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config211.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config211, source, orbit, processor))
    
    config212 = Config()
    config212.Algorithm_CSF_background_geometry = Geometry.CENTER
    config212.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config212.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config212, source, orbit, processor))

    config221 = Config() 
    config221.Algorithm_CSF_background_geometry = Geometry.CENTER
    config221.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config221.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config221, source, orbit, processor))

    config222 = Config()
    config222.Algorithm_CSF_background_geometry = Geometry.CENTER
    config222.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config222.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config222, source, orbit, processor))

    config231 = Config() 
    config231.Algorithm_CSF_background_geometry = Geometry.CENTER
    config231.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config231.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config231, source, orbit, processor))
    
    config232 = Config()
    config232.Algorithm_CSF_background_geometry = Geometry.CENTER
    config232.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config232.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config232, source, orbit, processor))

    config311 = Config() 
    config311.Algorithm_CSF_background_geometry = Geometry.INTERSECTS
    config311.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config311.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config311, source, orbit, processor))
    
    config312 = Config()
    config312.Algorithm_CSF_background_geometry = Geometry.INTERSECTS
    config312.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config312.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config312, source, orbit, processor))

    config321 = Config() 
    config321.Algorithm_CSF_background_geometry = Geometry.INTERSECTS
    config321.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config321.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config321, source, orbit, processor))

    config322 = Config()
    config322.Algorithm_CSF_background_geometry = Geometry.INTERSECTS
    config322.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config322.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config322, source, orbit, processor))

    config331 = Config() 
    config331.Algorithm_CSF_background_geometry = Geometry.INTERSECTS
    config331.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config331.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config331, source, orbit, processor))
    
    config332 = Config()
    config332.Algorithm_CSF_background_geometry = Geometry.INTERSECTS
    config332.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config332.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config332, source, orbit, processor))
    
    csv = []
    csv.append("Background geometry, Downwind box geometry, Downwind box mask, q (t/hr)")
    for m in models:
        if not path.exists(m):
            message = "Model {} is missing. Please run models first.".format(m)
            print(message)
            if not logger is None: logger.error(message)
            return

        with open(m, "rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()
        
        v = sad.q_valid_hour
        csv.append("{}, {}, {}, {}".format(sad.config.Algorithm_CSF_background_geometry.name, sad.config.Algorithm_CSF_downwindbox_geometry.name, sad.config.Algorithm_CSF_downwindbox_mask.name, v))

    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_comparison_downindbox_length(source:str, orbit: str, processor:str, fileout:str) -> None:
    ''' 
    Generate a csv file and print out comparison of emission reates for different configurations
    Note this calculation dows not depend on Geometry of background
    '''
    models = []
    
    config111 = Config() 
    config111.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config111.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config111.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config111, source, orbit, processor))
    
    config112 = Config()
    config112.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config112.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config112.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config112, source, orbit, processor))

    config121 = Config() 
    config121.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config121.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config121.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config121, source, orbit, processor))

    config122 = Config()
    config122.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config122.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config122.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config122, source, orbit, processor))

    config131 = Config() 
    config131.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config131.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config131.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config131, source, orbit, processor))
    
    config132 = Config()
    config132.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config132.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config132.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config132, source, orbit, processor))
    
    csv = []
    csv.append("Downwind box geometry, Downwind box mask, Length (km)")
    for m in models:
        if not path.exists(m):
            message = "Model {} is missing. Please run models first.".format(m)
            print(message)
            if not logger is None: logger.error(message)
            return

        with open(m, "rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()
        
        v = sad.downwind_box_length / 1000
        csv.append("{}, {}, {}".format(sad.config.Algorithm_CSF_downwindbox_geometry.name, sad.config.Algorithm_CSF_downwindbox_mask.name, v))

    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_comparison_downindbox_rotation(source:str, orbit: str, processor:str, fileout:str) -> None:
    ''' 
    Generate a csv file and print out comparison of emission reates for different configurations
    Note this calculation dows not depend on Geometry of background
    '''
    models = []
    
    config111 = Config() 
    config111.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config111.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config111.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config111, source, orbit, processor))
    
    config112 = Config()
    config112.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config112.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config112.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config112, source, orbit, processor))

    config121 = Config() 
    config121.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config121.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config121.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config121, source, orbit, processor))

    config122 = Config()
    config122.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config122.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config122.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config122, source, orbit, processor))

    config131 = Config() 
    config131.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config131.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config131.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config131, source, orbit, processor))
    
    config132 = Config()
    config132.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config132.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config132.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config132, source, orbit, processor))
    
    csv = []
    csv.append("Downwind box geometry, Downwind box mask, Rotation (deg)")
    for m in models:
        if not path.exists(m):
            message = "Model {} is missing. Please run models first.".format(m)
            print(message)
            if not logger is None: logger.error(message)
            return

        with open(m, "rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()
        
        v = sad.downwind_box_azimuth_adjusted - sad.downwind_box_azimuth_initial
        csv.append("{}, {}, {}".format(sad.config.Algorithm_CSF_downwindbox_geometry.name, sad.config.Algorithm_CSF_downwindbox_mask.name, v))

    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_comparison_downindbox_width(source:str, orbit: str, processor:str, fileout: str) -> None:
    ''' 
    Generate a csv file and print out comparison of emission rates for different configurations
    '''
    models = []
    
    config111 = Config() 
    config111.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config111.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config111.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config111, source, orbit, processor))
    
    config112 = Config()
    config112.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config112.Algorithm_CSF_downwindbox_geometry = Geometry.CONTAINS
    config112.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config112, source, orbit, processor))

    config121 = Config() 
    config121.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config121.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config121.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config121, source, orbit, processor))

    config122 = Config()
    config122.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config122.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    config122.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config122, source, orbit, processor))

    config131 = Config() 
    config131.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config131.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config131.Algorithm_CSF_downwindbox_mask = Mask.NONE
    models.append(Algorithm_CSF.get_picklename(config131, source, orbit, processor))
    
    config132 = Config()
    config132.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    config132.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    config132.Algorithm_CSF_downwindbox_mask = Mask.NEGATIVE
    models.append(Algorithm_CSF.get_picklename(config132, source, orbit, processor))

    csv = []
    csv.append("Downwind box geometry, Downwind box mask, Width (km)")
    for m in models:
        if not path.exists(m):
            message = "Model {} is missing. Please run models first.".format(m)
            print(message)
            if not logger is None: logger.error(message)
            return

        with open(m, "rb") as f:
            sad: Algorithm_CSF = load(f)
            f.close()
        
        v = sad.downwind_box_half_width * 2 / 1000
        csv.append("{}, {}, {}".format(sad.config.Algorithm_CSF_downwindbox_geometry.name, sad.config.Algorithm_CSF_downwindbox_mask.name, v))

    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_downindbox_rotation(orbit: str, processor: str, fileout: str) -> None:
    ''' 
    Generate a csv file of rotations of average enhancement for default configurations
    '''
    config = Config()
    m = Algorithm_CSF.get_picklename(config, "HailCreek", orbit, processor)
    
    csv = []
    csv.append("Rotation (deg), Azimuth (deg), Average Enhancement (ppb)")

    if not path.exists(m):
        message = "Model {} is missing. Please run models first.".format(m)
        print(message)
        if not logger is None: logger.error(message)
        return

    with open(m, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()
    
    arr = sad.rotations.copy()
    arr.sort(key=lambda x: x.rotation)
    
    for r in arr:
        if r.averageenhancement is None: csv.append("{}, {}, NaN".format(r.rotation, r.azimuth))
        else: csv.append("{}, {}, {}".format(r.rotation, r.azimuth, r.averageenhancement))

    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_downindbox_enhancementlength(orbit: str, processor: str, fileout:str) -> None:
    ''' 
    Generate a csv file and print out comparison of emission reates for different configurations
    Note this calculation dows not depend on Geometry of background
    '''
    config = Config()
    m = Algorithm_CSF.get_picklename(config, "HailCreek", orbit, processor)
    
    csv = []
    csv.append("Length (m), Total Enhancement (ppb), Increment")

    if not path.exists(m):
        message = "Model {} is missing. Please run models first.".format(m)
        print(message)
        if not logger is None: logger.error(message)
        return

    with open(m, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()
    
    arr = sad.enhancement_length.copy()
    
    for e in arr:
        if e.increment is None: csv.append("{}, {}, NaN".format(e.length, e.enhancement))
        else: csv.append("{}, {}, {}".format(e.length, e.enhancement, e.increment))

    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_downindbox_enhancementwidth(orbit: str, processor: str, fileout: str) -> None:
    ''' 
    Generate a csv file and print out comparison of emission reates for different configurations
    Note this calculation dows not depend on Geometry of background
    '''    
    m = Algorithm_CSF.get_picklename(Config(), "HailCreek", orbit, processor)
    
    csv = []
    csv.append("Half-width (m), Total Enhancement (ppb), Increment")

    if not path.exists(m):
        message = "Model {} is missing. Please run models first.".format(m)
        print(message)
        if not logger is None: logger.error(message)
        return

    with open(m, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()
    
    arr = sad.enhancement_width.copy()
    
    for e in arr:
        if e.increment is None: csv.append("{}, {}, NaN".format(e.half_width, e.enhancement))
        else: csv.append("{}, {}, {}".format(e.half_width, e.enhancement, e.increment))

    __utility_table_output_csv(fileout, csv)

    return

def _list_csf_transects(source:str, orbit: str, processor: str, fileout: str) -> None:
    '''
    List transects for current configuration of run for orbit
    '''
    config = Config()    
    filename = Algorithm_CSF.get_picklename(config, source, orbit, processor)
    if not path.exists(filename):
        print("File {} does not exist. Terminating".format(filename))
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    csv = ["id, valid, enhancement (kg/m), plume width (m), positive pixel count, q (t/hr)"]
    i = 4
    rotation = sad.downwind_box_azimuth_adjusted - sad.downwind_box_azimuth_initial
    cosadjstment = cos(rotation/180 * pi)
    windspeed = sqrt(sad.u * sad.u + sad.v * sad.v)
    conversion = 3600.0 / 1000.0 # conversion kg/s -> t/hour
    z = windspeed * cosadjstment * conversion
    for t in sad.transects:
        csv.append("{}, {}, {}, {}, {}, {}".format(i, t.isvalid(), t.enhkg, t.len_values_positive, t.pix_values_positive, t.enhkg * z))
        i += 1
    
    __utility_table_output_csv(fileout, csv)

    return

def _list_ef(fileout: str) -> None:
    __utility_table_output_csv(fileout, csv_ef)
    return

def _list_ef_IPCC(fileout: str) -> None:
    __utility_table_output_csv(fileout, csv_ef_ipcc)
    return

def _list_ime(source:str, orbit: str, processor: str, fileout: str) -> None:
    '''
    List transects for current configuration of run for orbit
    '''
    config = Config()    
    filename = Algorithm_CSF.get_picklename(config, source, orbit, processor)
    if not path.exists(filename):
        print("File {} does not exist. Terminating".format(filename))
        return
    
    with open(filename, "rb") as f:
        sad: Algorithm_CSF = load(f)
        f.close()

    csv = ["id, Box length (km), Time (h), Q (t/hr), Q (Gg/year)"]
    i = 3
    for t in sad.transects:
        i += 1
        if not t.isvalid(): continue
        ime = Algorithm_IME_factory(sad, sad.background, i)  
        csv.append("{}, {}, {}, {}, {}".format(ime.transectid, ime.length * 0.001, ime.residence_hours, ime.q_valid_hour, ime.q_valid_year))

    __utility_table_output_csv(fileout, csv)

    return

def _list_tropomi_correlation(orbit: str, processor:str, fileout:str) -> None:
    '''
    List correlation tables for domain and orbit

    Parameters
    ----------
    orbit: str
        Orbit for analysis
    '''
    csv=[]
    configdomain = Config()
    configdomain.Algorithm_CSF_Filter_R_domain = True
    configorbit = Config()
    configorbit.Algorithm_CSF_Filter_R_domain = False

    filterdomain = TROPOMI_Filter("r", configdomain, datetime(2018, 4, 30), datetime(2020, 1, 1))
    filterall = TROPOMI_Filter("r", configorbit, datetime(2018, 4, 30), datetime(2020, 1, 1))
    
    tropomi = TROPOMI_for_orbit(orbit, processor)
    r_CH4_albedo_SWIR_domain, r_CH4_albedo_NIR_domain, r_CH4_areosol_SWIR_domain, r_CH4_areosol_NIR_domain = filterdomain.calculate_correlation(tropomi)
    r_CH4_albedo_SWIR_all, r_CH4_albedo_NIR_all, r_CH4_areosol_SWIR_all, r_CH4_areosol_NIR_all = filterall.calculate_correlation(tropomi)
    
    csv.append("Orbit {}, Domain, Orbit".format(orbit))
    csv.append("R albedo SWIR, {}, {}". format(r_CH4_albedo_SWIR_domain, r_CH4_albedo_SWIR_all))
    csv.append("R albedo NIR, {}, {}". format(r_CH4_albedo_NIR_domain, r_CH4_albedo_NIR_all))
    csv.append("R aerosol SWIR, {}, {}". format(r_CH4_areosol_SWIR_domain, r_CH4_areosol_SWIR_all))
    csv.append("R aerosol NIR, {}, {}". format(r_CH4_areosol_NIR_domain, r_CH4_areosol_NIR_all))
    
    __utility_table_output_csv(fileout, csv)

    return

def _handler_chart_activity(args):
    _chart_HailCreek_Activity(None)

def _handler_chart_S2(args):
    if logger is None: Config.create_log("Paper_Run.log")
    _chart_sadavarte_figure2(args.type, None)

def _handler_chart_srtm(args):
    config = Config()
    _chart_SRTM(config.Algorithm_CSF_domain_ymax, config.Algorithm_CSF_domain_xmin, config.Algorithm_CSF_domain_ymin, config.Algorithm_CSF_domain_xmax, None)

def _handler_chart_synthetic(args):
    _chart_synthetic_transect(None)

def _handler_clean(args):
    '''
    Clean all outputs for algorithm
    '''
    Config.create_log("Paper_Clean.log")

    top = Config.Paper_folder
    if not path.exists(top):
        if not logger is None: logger.error("No folder {}".format(top))
        exit()

    papers = listdir(top)
    for p in papers:
        po = path.join(top,p)
        if path.isfile(po) and (p.startswith("Figure_") or p.startswith("Table_") or p.startswith("TROPOMI_Filter_")): 
            remove(po)
            print("Removed: {}".format(po))
            if not logger is None: logger.info("Removed: {}".format(po))

def _handler_list_csf_orbit(args):
    '''
    Handler for listing options
    '''

    if args.type == "background":
        _list_csf_comparison_background(args.source, args.orbit, None)
        return
    
    if args.type == "emission":
        _list_csf_comparison_emissionrate(args.source, args.orbit, None)
        return

    if args.type == "length":
        _list_csf_comparison_downindbox_length(args.source, args.orbit, None)
        _list_csf_downindbox_enhancementlength(args.orbit, None)
        return
    
    if args.type == "rotation":
        _list_csf_comparison_downindbox_rotation(args.source, args.orbit, None)
        _list_csf_downindbox_rotation(args.orbit, None)
        return
    
    if args.type == "transect":
        _list_csf_transects(args.source, args.orbit, None)
        return

    if args.type == "width":
        _list_csf_comparison_downindbox_width(args.source, args.orbit, None)
        _list_csf_downindbox_enhancementwidth(args.orbit, None)
        return
    
    print("There is no list for this selection source:{} type:{}".format(args.source, args.type))

def _handler_list_csf_summary(args):
    '''
    Handler for listing options
    '''  
    if args.type == "figure2_all":
        _list_csf_Figure2("figure2_all", None)
        return

    if args.type == "figure2_min_max":
        _list_csf_Figure2("figure2_min_max", None)
        return

    print("There is no list for this selection type:{}".format(args.type))

def _handler_list_era5(args):
    '''
    Handler for listing options
    '''
    if args.type == "blh":
        _list_ERA5_BLH("HailCreek", "09956")
        exit()
    if args.type == "sp":
        _list_ERA5_SP("HailCreek", "09956")
        exit()
    if args.type == "wind":
        _list_ERA5_Wind("HailCreek", "09956", None)
        exit()

    print("There is no list for this selection source:{} type:{}".format(args.source, args.type))

def _handler_list_ime(args):
    _list_ime(args.source, args.orbit, path.join(Config.Paper_folder, None))


def _handler_list_tm(args):
    '''
    Handler for listing options for Integrated Mass Enhancement method
    '''
    _list_csf_tm_comparison(args.source, args.orbit, args.minbackground, args.maxbackground, None)
    return

    # print("There is no list for this selection source:{} type:{}".format(args.source, args.type))

def _handler_list_tropomi(args):
    '''
    Handler for listing options
    '''
    if args.type == "correlations":
        _list_tropomi_correlation("03912", None)
        _list_tropomi_correlation("02990", None)
        return
    
    if args.type in ["contains", "qa", "r"]:
        _list_tropomi_filteredorbits(args.type,None)
        return
        
    print("There is no list for this selection source:{} type:{}".format(args.source, args.type))

def _handler_run(args):
    '''
    Run everything needed all charts and tables used in the paper
    '''
    Config.create_log("Paper_Run.log")
    source = create_source(args.source)

    # Run all calculations used by Paper
    # Run CSF for orbit 02990 default configuration. Note that this orbit is rejected by TROPOMI filters
    config = Config()
    orbits = []
    with open(TROPOMI_Filter.pickle, "rb") as f:
        filter: TROPOMI_Filter = load(f)
        orbits = filter.result.copy()
        f.close()
    ncdict = filter.find_highest_configured_procesor(config.TROPOMI_Processors)  
    processor2990 = TROPOMI.get_processor_version(ncdict[2990])

    filename = Algorithm_CSF.get_picklename(config, "HailCreek", "02990", processor2990)
    if not path.exists(filename):
        tropomi = TROPOMI_for_orbit("02990", processor2990)
        [_,_,dt,_] = tropomi.get_pixel_for_source(source.xy)
        grid_time = ERA5.calculate_ERA_GridTime(dt)

        path_ERA5_SL = ERA5_SingleLevel.get_filepath("CSF", args.source, "02990")
        if (not(path.exists(path_ERA5_SL))):
            message = "ERA5 data {} is missing. Please use ERA5 program to download first.".format(path_ERA5_SL)
            if not logger is None: logger.error(message)
            print(message)
            return

        era5_sl = ERA5_SingleLevel.from_netCDF(path_ERA5_SL)
        
        path_ERA5_PL = ERA5_PressureLevels.get_filepath("CSF",source.case_name, "02990")
        if not(path.exists(path_ERA5_PL)) :
            message = "ERA5 data {} is missing. Please use ERA5 program to download first.".format(path_ERA5_PL)
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

        sad = Algorithm_CSF.run(config, source, "02990", processor2990, wind)

    # Run 09956 for different gemetries and masking configurations
    processor9956 = TROPOMI.get_processor_version(ncdict[9956])
    tropomi = TROPOMI_for_orbit("09956", processor9956)
    [_,_,dt,_] = tropomi.get_pixel_for_source(source.xy)
    grid_time = ERA5.calculate_ERA_GridTime(dt)

    path_ERA5_SL = ERA5_SingleLevel.get_filepath("CSF", args.source, "09956")
    if (not(path.exists(path_ERA5_SL))):
        message = "ERA5 data {} is missing. Please use ERA5 program to download first.".format(path_ERA5_SL)
        if not logger is None: logger.error(message)
        print(message)
        return

    era5_sl = ERA5_SingleLevel.from_netCDF(path_ERA5_SL)
    
    path_ERA5_PL = ERA5_PressureLevels.get_filepath("CSF",source.case_name, "09956")
    if not(path.exists(path_ERA5_PL)) :
        message = "ERA5 data {} is missing. Please use ERA5 program to download first.".format(path_ERA5_PL)
        if not logger is None: logger.error(message)
        print(message)
        return

    era5_p = ERA5_PressureLevels.from_netCDF(path_ERA5_PL)

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

    for bg in Geometry:
        for dg in Geometry:
            for m in Mask:
                config = Config()
                config.Algorithm_CSF_background_geometry = bg
                config.Algorithm_CSF_downwindbox_geometry = dg
                config.Algorithm_CSF_downwindbox_mask = m
                filename = Algorithm_CSF.get_picklename(config, source.case_name, "09956", processor9956)
                if not path.exists(filename):
                    _ = Algorithm_CSF.run(config, source, "09956", processor9956, wind)

    # Run all orbits passing R filter with default configuration
    for o in orbits:
        if not (o in ncdict): continue
        processor = TROPOMI.get_processor_version(ncdict[o])
        strorbit = str(o).zfill(5)
        filename = Algorithm_CSF.get_picklename(config, source.case_name, strorbit, processor)
        if path.exists(filename): continue

        tropomi = TROPOMI_for_orbit(strorbit, processor)
        [_,_,dt,_] = tropomi.get_pixel_for_source(source.xy)
        grid_time = ERA5.calculate_ERA_GridTime(dt)
        
        path_ERA5_SL = ERA5_SingleLevel.get_filepath("CSF", source.case_name, strorbit)
        if (not(path.exists(path_ERA5_SL))):
            message = "Required file: {} does not exist. Please download this file first using ERA5 program.".format(path_ERA5_SL)
            if not logger is None: logger.error(message)
            print(message)
            continue

        era5_sl = ERA5_SingleLevel.from_netCDF(path_ERA5_SL)
        
        path_ERA5_PL = ERA5_PressureLevels.get_filepath("CSF",source.case_name, strorbit)
        if not(path.exists(path_ERA5_PL)) :
            message = "Required file: {} does not exist. Please download this file first using ERA5 program.".format(path_ERA5_PL)
            if not logger is None: logger.error(message)
            print(message)
            continue

        era5_p = ERA5_PressureLevels.from_netCDF(path_ERA5_PL)

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

        config = Config()
        _ = Algorithm_CSF.run(config, source, strorbit, processor, wind)

    orbit = "09956"
    model = "GFSQ"
    emission_start_date = datetime(2019, 9, 14, 4, tzinfo=timezone.utc)
    emission_end_date = datetime(2019, 9, 15, 3, tzinfo=timezone.utc)
    chart_datetime = datetime(2019, 9, 15, 4, tzinfo=timezone.utc)
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Trajectory_09956_2019091404_2019091503_2019091504.png")

    _chart_HYSPLIT_HailCreek_SingleTrajectory("HailCreek", orbit, processor9956, model, emission_start_date, emission_end_date, chart_datetime, fileout)
    _chart_SRTM(-20.5, 146.0, -24.0, 150.0, path.join(Config.Paper_folder, "Figure_SRTM_Topography_HailCreek_Large.png"))
    _chart_SRTM(-21.0, 147.85, -22.0, 148.9, path.join(Config.Paper_folder, "Figure_SRTM_Topography_HailCreek_Detail.png"))
    _chart_AWS_windrose(path.join(Config.Paper_folder, "Figure_MoranbahAWS_windrose.png"))
    _chart_HailCreek_Activity(path.join(Config.Paper_folder, "Figure_CoalProduction_HailCreek_FinancialYear.png"))
    _chart_synthetic_geometry(path.join(Config.Paper_folder,"Figure_Synthetic_Geometry.png"))
    _chart_synthetic_masking(path.join(Config.Paper_folder,"Figure_Synthetic_Mask.png"))
    _chart_synthetic_transect(path.join(Config.Paper_folder,"Figure_Synthetic_Background.png"))
    _chart_CSF_background_box_combined("09956", processor9956, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09956_BackgroundComparison.png"))
    _chart_CSF_elements("HailCreek", "09956", processor9956, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09956_Elements.png"))

    c_14 = Config()
    c_14.Algorithm_CSF_background_geometry = Geometry.CONTAINS
    c_14.Algorithm_CSF_downwindbox_geometry = Geometry.INTERSECTS
    c_14.Algorithm_CSF_downwindbox_mask = Mask.NONE
    _chart_CSF_positive_config("HailCreek", "09956", processor9956, c_14, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09956_CONTAINS_INTERSECT_NONE_Positive.png"))

    c_15 = Config()
    c_15.Algorithm_CSF_background_geometry = Geometry.INTERSECTS
    c_15.Algorithm_CSF_downwindbox_geometry = Geometry.CENTER
    c_15.Algorithm_CSF_downwindbox_mask = Mask.NONE
    _chart_CSF_positive_config("HailCreek", "09956", processor9956, c_15, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09956_INTERSECT_CENTER_NONE_Positive.png"))

    TROPOMI_for_orbit("02990", processor2990).chart_size(source, Geometry.INTERSECTS, 1.0, path.join(Config.Paper_folder, "Figure_TROPOMI_02990.png"))
    _chart_CSF_valid("HailCreek", "02990", processor2990, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_02990_Valid.png"))

    processor6906 = TROPOMI.get_processor_version(ncdict[6906])
    TROPOMI_for_orbit("06906", processor6906).chart_size(source, Geometry.INTERSECTS, 1.0, path.join(Config.Paper_folder, "Figure_TROPOMI_06906.png"))
    _chart_CSF_valid("HailCreek", "06906", processor6906, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_06906_Valid.png"))
    _chart_CSF_positive("HailCreek", "06906", processor6906, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_06906_Positive.png"))

    processor11020 = TROPOMI.get_processor_version(ncdict[11020])
    TROPOMI_for_orbit("11020", processor11020).chart_size(source, Geometry.INTERSECTS, 1.0, path.join(Config.Paper_folder, "Figure_TROPOMI_11020.png"))
    _chart_CSF_valid("HailCreek", "11020", processor11020, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_11020_Valid.png"))
    
    processor9488 = TROPOMI.get_processor_version(ncdict[9488])
    TROPOMI_for_orbit("09488", processor9488).chart_size(source, Geometry.INTERSECTS, 1.0, path.join(Config.Paper_folder, "Figure_TROPOMI_09488.png"))
    _chart_CSF_valid("HailCreek", "09488", processor9488, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09488_Valid.png"))

    _chart_sadavarte_figure2("figure2_all", path.join(Config.Paper_folder, "Figure_CSF_Figure2_All.png"))
    _chart_sadavarte_figure2("figure2_min_max", path.join(Config.Paper_folder,"Figure_CSF_Figure2_Filtered_Min_Max.png"))

    ################# Case study for orbit 09956

    chart_MSLP("HailCreek", "09956", processor9956, path.join(Config.Paper_folder, "Figure_MSLP_09956.png"))
    TROPOMI_for_orbit("09956", processor9956).chart_size(source, Geometry.INTERSECTS, 1.0, path.join(Config.Paper_folder, "Figure_TROPOMI_09956.png"))
    _chart_CSF_valid("HailCreek", "09956", processor9956, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09956_Valid.png"))
    _chart_CSF_positive("HailCreek", "09956",  processor9956, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09956_Positive.png"))
    _chart_CSF_transects("HailCreek", "09956", processor9956, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09956_Transect.png"))
    
    filter_start_start = datetime(2019, 9, 14, 4, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 9, 15, 3, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 9, 15, 4, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 9, 15, 4, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 9, 15, 4, tzinfo=timezone.utc)
    orbit = "09956"
    model =  "GFSQ"
    direction = "F"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_09956_GFSQ_F_2019091404_2019091503_2019091504_2019091504_2019091504.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor9956, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)

    filter_start_start = datetime(2019, 9, 14, 9, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 9, 15, 3, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 9, 15, 4, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 9, 15, 4, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 9, 15, 4, tzinfo=timezone.utc) 
    orbit = "09956"
    model =  "GFSQ"
    direction = "F"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_09956_GFSQ_F_2019091409_2019091503_2019091504_2019091504_2019091504.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor9956, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)

    processor9942 = TROPOMI.get_processor_version(ncdict[9942])
    TROPOMI_for_orbit("09942", processor9942).chart_size(source, Geometry.INTERSECTS, 1.0, path.join(Config.Paper_folder, "Figure_TROPOMI_09942.png"))

    filter_start_start = datetime(2019, 9, 14, 4, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 9, 15, 3, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 9, 14, 5, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 9, 15, 3, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 9, 14, 4, tzinfo=timezone.utc)
    orbit = "09942"
    model =  "GFSQ"
    direction = "B"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_09942_GFSQ_B_2019091404_2019091503_2019091405_2019091503_2019091404.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor9942, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)

    filter_start_start = datetime(2019, 9, 14, 4, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 9, 14, 4, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 9, 14, 5, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 9, 14, 9, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 9, 14, 4, tzinfo=timezone.utc)
    orbit = "09942"
    model =  "GFSQ"
    direction = "B"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_09942_GFSQ_B_2019091404_2019091404_2019091405_2019091409_2019091404.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor9942, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)

    _chart_TM("HailCreek", "09956", processor9956, "GFSQ", "F", 1812, path.join(Config.Paper_folder,"Figure_TM_09956_1812.png"))
    _chart_TM("HailCreek", "09956", processor9956, "GFSQ", "F", 1815, path.join(Config.Paper_folder,"Figure_TM_09956_1815.png"))
    _chart_TM("HailCreek", "09956", processor9956, "GFSQ", "F", 1820, path.join(Config.Paper_folder,"Figure_TM_09956_1820.png"))

    _chart_IME("HailCreek", "09956", processor9956, 1812, 4, path.join(Config.Paper_folder,"Figure_IME_09956_1812_4.png"))
    _chart_IME("HailCreek", "09956", processor9956, 1812, 8, path.join(Config.Paper_folder,"Figure_IME_09956_1812_8.png"))
    _chart_IME("HailCreek", "09956", processor9956, 1812, 12, path.join(Config.Paper_folder,"Figure_IME_09956_1812_12.png"))
    
    _chart_IME("HailCreek", "09956", processor9956, 1815, 4, path.join(Config.Paper_folder,"Figure_IME_09956_1815_4.png"))
    _chart_IME("HailCreek", "09956", processor9956, 1815, 8, path.join(Config.Paper_folder,"Figure_IME_09956_1815_8.png"))
    _chart_IME("HailCreek", "09956", processor9956, 1815, 12, path.join(Config.Paper_folder,"Figure_IME_09956_1815_12.png"))
    
    _chart_IME("HailCreek", "09956", processor9956, 1820, 4, path.join(Config.Paper_folder,"Figure_IME_09956_1820_4.png"))
    _chart_IME("HailCreek", "09956", processor9956, 1820, 8, path.join(Config.Paper_folder,"Figure_IME_09956_1820_8.png"))
    _chart_IME("HailCreek", "09956", processor9956, 1820, 12, path.join(Config.Paper_folder,"Figure_IME_09956_1820_12.png"))
    
    ################# Case study for orbit 11332
    processor11332 = TROPOMI.get_processor_version(ncdict[11332])
    chart_MSLP("HailCreek", "11332", processor11332, path.join(Config.Paper_folder, "Figure_MSLP_11332.png"))
    TROPOMI_for_orbit("11332", processor11332).chart_size(source, Geometry.INTERSECTS, 1.0, path.join(Config.Paper_folder, "Figure_TROPOMI_11332.png"))
    _chart_CSF_valid("HailCreek", "11332", processor11332, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_11332_Valid.png"))
    _chart_CSF_positive("HailCreek", "11332", processor11332, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_11332_Positive.png"))
    _chart_CSF_transects("HailCreek", "11332", processor11332, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_11332_Transect.png"))

    filter_start_start = datetime(2019, 12, 20, 4, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 12, 21, 3, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 12, 21, 4, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 12, 21, 4, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 12, 21, 4, tzinfo=timezone.utc)    
    orbit = "11332"
    model =  "GFSQ"
    direction = "F"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_11332_GFSQ_F_2019122004_2019122103_2019122104_2019122104_2019122104.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor11332, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)

    
    processor11318 = TROPOMI.get_processor_version(ncdict[11318])
    filter_start_start = datetime(2019, 12, 20, 4, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 12, 20, 4, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 12, 20, 4, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 12, 21, 3, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 12, 20, 4, tzinfo=timezone.utc)    
    orbit = "11318"
    model =  "GFSQ"
    direction = "B"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_11318_GFSQ_B_2019122004_2019122004_2019122004_2019122103_2019122004.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor11318, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)
    
    filter_start_start = datetime(2019, 12, 20, 4, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 12, 20, 4, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 12, 20, 5, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 12, 20, 9, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 12, 20, 4, tzinfo=timezone.utc)    
    orbit = "11318"
    model =  "GFSQ"
    direction = "B"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_11318_GFSQ_B_2019122004_2019122004_2019122005_2019122009_2019122004.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor11318, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)

    _chart_TM("HailCreek", "11332", processor11332, "GFSQ", "F", 1809, path.join(Config.Paper_folder,"Figure_TM_11332_1809.png"))
    _chart_TM("HailCreek", "11332", processor11332, "GFSQ", "F", 1812, path.join(Config.Paper_folder,"Figure_TM_11332_1812.png"))
    _chart_TM("HailCreek", "11332", processor11332, "GFSQ", "F", 1815, path.join(Config.Paper_folder,"Figure_TM_11332_1815.png"))

    _chart_IME("HailCreek", "11332", processor11332, 1809, 4, path.join(Config.Paper_folder,"Figure_IME_11332_1809_4.png"))
    _chart_IME("HailCreek", "11332", processor11332, 1809, 8, path.join(Config.Paper_folder,"Figure_IME_11332_1809_8.png"))
    _chart_IME("HailCreek", "11332", processor11332, 1809, 12, path.join(Config.Paper_folder,"Figure_IME_11332_1809_12.png"))
    _chart_IME("HailCreek", "11332", processor11332, 1809, 15, path.join(Config.Paper_folder,"Figure_IME_11332_1809_15.png"))

    _chart_IME("HailCreek", "11332", processor11332, 1812, 4, path.join(Config.Paper_folder,"Figure_IME_11332_1812_4.png"))
    _chart_IME("HailCreek", "11332", processor11332, 1812, 8, path.join(Config.Paper_folder,"Figure_IME_11332_1812_8.png"))
    _chart_IME("HailCreek", "11332", processor11332, 1812, 12, path.join(Config.Paper_folder,"Figure_IME_11332_1812_12.png"))
    _chart_IME("HailCreek", "11332", processor11332, 1812, 15, path.join(Config.Paper_folder,"Figure_IME_11332_1812_15.png"))

    _chart_IME("HailCreek", "11332", processor11332, 1815, 4, path.join(Config.Paper_folder,"Figure_IME_11332_1815_4.png"))
    _chart_IME("HailCreek", "11332", processor11332, 1815, 8, path.join(Config.Paper_folder,"Figure_IME_11332_1815_8.png"))
    _chart_IME("HailCreek", "11332", processor11332, 1815, 12, path.join(Config.Paper_folder,"Figure_IME_11332_1815_12.png"))
    _chart_IME("HailCreek", "11332", processor11332, 1815, 15, path.join(Config.Paper_folder,"Figure_IME_11332_1815_15.png"))

    ################# Case study for orbit 09445
    processor9445 = TROPOMI.get_processor_version(ncdict[9445])
    chart_MSLP("HailCreek", "09445", processor9445, path.join(Config.Paper_folder, "Figure_MSLP_09445.png"))
    TROPOMI_for_orbit("09445", processor9445).chart_size(source, Geometry.INTERSECTS, 1.0, path.join(Config.Paper_folder, "Figure_TROPOMI_09445.png"))
    _chart_CSF_valid("HailCreek", "09445", processor9445, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09445_Valid.png"))
    _chart_CSF_positive("HailCreek", "09445", processor9445, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09445_Positive.png"))
    _chart_CSF_transects("HailCreek", "09445", processor9445, path.join(Config.Paper_folder, "Figure_CSF_HailCreek_09445_Transect.png"))

    filter_start_start = datetime(2019, 8, 10, 0, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 8, 10, 3, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 8, 10, 3, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 8, 10, 3, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 8, 10, 3, tzinfo=timezone.utc)    
    orbit = "09445"
    model =  "GFSQ"
    direction = "F"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_09445_GFSQ_F_2019081000_2019081003_2019081003_2019081003_2019081003.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor9445, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)
    
    processor9431 = TROPOMI.get_processor_version(ncdict[9431])
    filter_start_start = datetime(2019, 8, 9, 4, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 8, 9, 4, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 8, 9, 4, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 8, 10, 3, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 8, 9, 4, tzinfo=timezone.utc)    
    orbit = "09431"
    model =  "GFSQ"
    direction = "B"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_09431_GFSQ_B_2019080904_2019080904_2019080904_2019081003_2019080904.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor9431, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)

    filter_start_start = datetime(2019, 8, 9, 4, tzinfo=timezone.utc)
    filter_start_end = datetime(2019, 8, 9, 4, tzinfo=timezone.utc)
    filter_end_start = datetime(2019, 8, 9, 23, tzinfo=timezone.utc)
    filter_end_end = datetime(2019, 8, 10, 3, tzinfo=timezone.utc)
    slice_datetime = datetime(2019, 8, 9, 4, tzinfo=timezone.utc)
    orbit = "09431"
    model =  "GFSQ"
    direction = "B"
    fileout = path.join(Config.Paper_folder, "Figure_HYSPLIT_Slice_09431_GFSQ_B_2019080904_2019080904_2019080923_2019081003_2019080904.png")
    _chart_HYSPLIT_time_slice("HailCreek", orbit, processor9431, model, direction, filter_start_start, filter_start_end, filter_end_start, filter_end_end, slice_datetime, fileout)

    # Generate tables

    _list_tropomi_filteredorbits("contains", path.join(Config.Paper_folder, "Table_TROPOMIFilter_contains.csv"))
    _list_tropomi_filteredorbits("qa", path.join(Config.Paper_folder, "Table_TROPOMIFilter_qa.csv"))
    _list_tropomi_filteredorbits("r", path.join(Config.Paper_folder, "Table_TROPOMIFilter_r.csv"))

    _list_csf_Figure2("figure2_all", path.join(Config.Paper_folder, "Table_CSF_HailCreek_Figure2_All.csv"))
    _list_csf_Figure2("figure2_min_max", path.join(Config.Paper_folder, "Table_CSF_HailCreek_Figure2_Filtered_Min_Max.csv"))

    processor3912 = TROPOMI.get_processor_version(ncdict[3912])
    _list_ef(path.join(Config.Paper_folder, "Table_IEF.csv"))
    _list_ef_IPCC(path.join(Config.Paper_folder, "Table_IEF_IPCC.csv"))
    _list_HailCreek(path.join(Config.Paper_folder, "Table_HailCreek.csv"))
    _list_activity_HailCreek_CY(path.join(Config.Paper_folder, "Table_HailCreek_Production_2025-calendar-year-coal-production-statistics.csv"))
    _list_activity_HailCreek_FY(path.join(Config.Paper_folder, "Table_HailCreek_Production_fy2010_fy2025.csv"))
    _list_HailCreek_emissions(path.join(Config.Paper_folder, "Table_HailCreek_emissions.csv"))
    _list_tropomi_correlation("03912", processor3912, path.join(Config.Paper_folder, "Table_TROPOMI_03912_correlation.csv"))
    _list_tropomi_correlation("02990", processor2990, path.join(Config.Paper_folder, "Table_TROPOMI_02990_correlation.csv"))
    _list_ERA5_Wind("HailCreek", "09956",  processor9956, path.join(Config.Paper_folder, "Table_ERA5_09956_Wind.csv"))

    date_start = datetime(2019,9,14,4,0, tzinfo=timezone.utc)
    date_end = datetime(2019,9,15,4,0, tzinfo=timezone.utc)
    _list_BOM_Moranbabah(date_start, date_end, path.join(Config.Paper_folder, "Table_BOM_Moranbah_2019091404_2019091504.csv"))
    _list_csf_downindbox_rotation("09956", processor9956, path.join(Config.Paper_folder, "Table_CSF_HailCreek_09956_Rotation.csv"))
    _list_csf_downindbox_enhancementlength("09956", processor9956, path.join(Config.Paper_folder, "Table_CSF_HailCreek_09956_EnhancementLength.csv"))
    _list_csf_downindbox_enhancementwidth("09956", processor9956, path.join(Config.Paper_folder, "Table_CSF_HailCreek_09956_EnhancementWidth.csv"))
    _list_csf_transects("HailCreek", "09956", processor9956, path.join(Config.Paper_folder, "Table_CSF_HailCreek_09956_Transect.csv"))

    # Table 20 Top 7 success orbits sorted ascending by number of pixels used for transect calculations
    # Table 21 Top 7 success orbits sorted descending by number of pixels used for transect calculations

    _list_csf_comparison_background("HailCreek","09956", processor9956, path.join(Config.Paper_folder, "Table_CSF_HailCreek_09956_Comparison_Background.csv")) 
    _list_csf_comparison_downindbox_rotation("HailCreek", "09956", processor9956, path.join(Config.Paper_folder, "Table_CSF_HailCreek_09956_Comparison_Rotation.csv"))
    _list_csf_comparison_downindbox_length("HailCreek", "09956", processor9956, path.join(Config.Paper_folder, "Table_CSF_HailCreek_09956_Comparison_Length.csv"))
    _list_csf_comparison_downindbox_width("HailCreek", "09956", processor9956, path.join(Config.Paper_folder, "Table_CSF_HailCreek_09956_Comparison_Width.csv"))
    _list_csf_comparison_emissionrate("HailCreek", "09956", processor9956, path.join(Config.Paper_folder, "Table_CSF_HailCreek_09956_Comparison_EmissionRate.csv"))
    _list_csf_tm_comparison("HailCreek", "09956", processor9956, 1800, 1825, path.join(Config.Paper_folder, "Table_Comparison_CSF_TM_HailCreek_09956.csv"))

    _list_ime("HailCreek", "09956",  processor9956, path.join(Config.Paper_folder, "Table_IME_HailCreek_09956.csv"))

    # Case study 11332
    date_start = datetime(2019,12,20,4,0, tzinfo=timezone.utc)
    date_end = datetime(2019,12,21,4,0, tzinfo=timezone.utc)
    _list_BOM_Moranbabah(date_start, date_end, path.join(Config.Paper_folder, "Table_BOM_Moranbah_2019122004_2019122104.csv"))
    _list_csf_tm_comparison("HailCreek", "11332", processor11332, 1800, 1825, path.join(Config.Paper_folder, "Table_Comparison_CSF_TM_HailCreek_11332.csv"))
    _list_ime("HailCreek", "11332",  processor11332, path.join(Config.Paper_folder, "Table_IME_HailCreek_11332.csv"))
    
    # Case study 09445
    date_start = datetime(2019,8,9,4,0, tzinfo=timezone.utc)
    date_end = datetime(2019,8,10,4,0, tzinfo=timezone.utc)
    _list_BOM_Moranbabah(date_start, date_end, path.join(Config.Paper_folder, "Table_BOM_Moranbah_2019080904_2019081004.csv"))

    # Table 27. Location of trajectory start points

if __name__ == '__main__':
    parser = ArgumentParser(prog="Paper")
    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_clean_algorithm = subparsers.add_parser("clean", help="Clean outputs for an algorithm.")
    parser_clean_algorithm.set_defaults(func=_handler_clean)

    parser_chart = subparsers.add_parser("chart", help="Chart.")
    sub_parser_chart = parser_chart.add_subparsers(help="subcommand help", required=True)

    parser_chart_activity= sub_parser_chart.add_parser("activity", help="Chart activity.")
    parser_chart_activity.set_defaults(func=_handler_chart_activity)

    parser_chart_S2 = sub_parser_chart.add_parser("S2", help="Chart S2 values.")
    parser_chart_S2.add_argument("type", type=str, choices=["figure2", "figure2_min_max"], help="List different data sects for paper")
    parser_chart_S2.set_defaults(func=_handler_chart_S2)

    parser_chart_SRTM = sub_parser_chart.add_parser("SRTM", help="Chart SRTM chart over domain.")
    parser_chart_SRTM.set_defaults(func=_handler_chart_srtm)

    parser_chart_synthetic = sub_parser_chart.add_parser("synthetic", help="Chart synthetic charts.")
    parser_chart_synthetic.set_defaults(func=_handler_chart_synthetic)

    # Parsers for listing tabular data used in paper
    parser_list = subparsers.add_parser("list", help="List commands.")
    subparsers_list = parser_list.add_subparsers(help="subcommand help", required=True)

    parser_list_csf = subparsers_list.add_parser("csf", help="List CSF data")
    subparsers_list_csf = parser_list_csf.add_subparsers(help="subcommand help", required=True)

    parser_list_csf_orbit = subparsers_list_csf.add_parser("orbit", help="List summary for all model runs")
    parser_list_csf_orbit.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_list_csf_orbit.add_argument("orbit", type=str, help="Orbit as 5 digit string e.g. 09956")
    parser_list_csf_orbit.add_argument("type", type=str, choices=["background", "emission", "tm", "length", "rotation", "transect", "width"], help="List different data sects for paper")
    parser_list_csf_orbit.set_defaults(func=_handler_list_csf_orbit)

    parser_list_csf_summary = subparsers_list_csf.add_parser("summary", help="List summary for all model runs")
    parser_list_csf_summary.add_argument("type", type=str, choices=["figure2_all","figure2_min_max"], help="List different data sects for paper")
    parser_list_csf_summary.set_defaults(func=_handler_list_csf_summary)

    parser_list_era5 = subparsers_list.add_parser("era5", help="List ERA5 data")
    parser_list_era5.add_argument("type", type=str, choices=["blh", "sp", "wind"], help="boundary layer height, surface pressure for domain 2019-09-14 00Z to 2019-09-15 04Z.")
    parser_list_era5.set_defaults(func=_handler_list_era5)

    parser_list_ime = subparsers_list.add_parser("ime", help="List IME data")
    parser_list_ime.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_list_ime.add_argument("orbit", type=str, help="Orbit as 5 digit string e.g. 09956")
    parser_list_ime.set_defaults(func=_handler_list_ime)

    parser_list_tm = subparsers_list.add_parser("tm", help="List TM data")
    parser_list_tm.add_argument("orbit", type=str, help="Orbit as 5 digit string e.g. 09956")
    parser_list_tm.add_argument("minbackground", type=float, help="Minimum background concentration inclusinve in ppb e.g. 1800")
    parser_list_tm.add_argument("maxbackground", type=float, help="Maximum background concentration exclusinve in ppb e.g. 1825")
    parser_list_tm.set_defaults(func=_handler_list_tm)

    parser_list_tropomi_filter = subparsers_list.add_parser("tropomi", help="List options relevane to TROPOMI")
    parser_list_tropomi_filter.add_argument("type", type=str, choices=["correlations", "contains", "qa", "r"], help="List different data sects for paper")
    parser_list_tropomi_filter.set_defaults(func=_handler_list_tropomi)

    # End Parsers for listing tabular data used in paper

    # Run all models and generate chart and tables used in this paper
    parser_run = subparsers.add_parser("run", help="Run document generation.")
    parser_run.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_run.set_defaults(func=_handler_run)

    args = parser.parse_args()    
    argumentsvalid = True

    if ("orbit" in args):
        if not(len(args.orbit) == 5) or not(args.orbit.isdigit()):
            print ("Orbit must be a five digit number: {}".format(args.orbit))
            argumentsvalid = False
    
    if argumentsvalid: args.func(args)