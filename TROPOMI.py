from argparse import ArgumentParser, Namespace
import cartopy.crs as ccrs
from BOM_AWS import create_aws
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime, timedelta, timezone
from Config import Config
from Geometry import Geometry
import logging
import matplotlib.pyplot as plt
from netCDF4 import Dataset, Variable
import numpy as np
from os import listdir, path
from shapely import box, MultiPolygon, Point, Polygon, prepare
from Source import Source, create_source

logger = logging.getLogger(__name__)

class TROPOMI:
    '''
    TROPOMI files downloaded from NCI THREDD server
    Please note:
    Until and including orbit 05832 from 2018-11-28 we are dealing with RPRO files
    Starting from orbit 05833 on 2018-11-28 we are dealing with OFFL files
    Possibly these should be separated into different objects
    '''
    # Dictionary of TROPOMI dimensions
    SP5CH4dimensions = {
        "scanline" : "/PRODUCT",
        "ground_pixel" : "/PRODUCT",
        "corner" : "/PRODUCT",
        "time" : "/PRODUCT",
        "layer" : "/PRODUCT",
        "level" : "/PRODUCT",
        "vertices" : "/METADATA/QA_STATISTICS",
        "XCH4_histogram_axis" : "/METADATA/QA_STATISTICS",
        "XCH4_pdf_axis(400)" : "/METADATA/QA_STATISTICS",
    }

    # Dictionary of TROPOMI variables
    SP5CH4variables = { 
        "scanline" : "/PRODUCT",
        "ground_pixel" : "/PRODUCT",
        "time" : "/PRODUCT",
        "corner" : "/PRODUCT",
        "layer" : "/PRODUCT",
        "level" : "/PRODUCT",
        "delta_time" : "/PRODUCT",
        "time_utc" : "/PRODUCT",
        "qa_value" : "/PRODUCT",
        "latitude" : "/PRODUCT",
        "longitude" : "/PRODUCT",
        "methane_mixing_ratio" : "/PRODUCT",
        "methane_mixing_ratio_precision" : "/PRODUCT",
        "methane_mixing_ratio_bias_corrected" : "/PRODUCT",

        "satellite_latitude" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "satellite_longitude" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "satellite_altitude" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "satellite_orbit_phase" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "solar_zenith_angle" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "solar_azimuth_angle" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "viewing_zenith_angle" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "viewing_azimuth_angle" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "latitude_bounds" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "longitude_bounds" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",
        "geolocation_flags" : "/PRODUCT/SUPPORT_DATA/GEOLOCATIONS",

        "processing_quality_flags" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "number_of_spectral_points_in_retrieval" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "number_of_spectral_points_in_retrieval_NIR" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "column_averaging_kernel" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "carbonmonoxide_total_column" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "carbonmonoxide_total_column_precision" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "water_total_column" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "water_total_column_precision" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "aerosol_size" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "aerosol_size_precision" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "aerosol_number_column" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "aerosol_number_column_precision" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "aerosol_mid_altitude" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "aerosol_mid_altitude_precision" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "surface_albedo_SWIR" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "surface_albedo_SWIR_precision" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "surface_albedo_NIR" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "surface_albedo_NIR_precision" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "aerosol_optical_thickness_SWIR" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "aerosol_optical_thickness_NIR" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "wavelength_calibration_offset_SWIR" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "wavelength_calibration_offset_NIR" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "chi_square" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "chi_square_SWIR" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "chi_square_NIR" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "degrees_of_freedom" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "degrees_of_freedom_methane" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "degrees_of_freedom_aerosol" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "number_of_iterations" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
        "fluorescence" : "/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",

        "altitude_levels" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "apparent_scene_pressure" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "apparent_scene_pressure_standard_deviation" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "cloud_fraction_VIIRS_NIR_IFOV" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "cloud_fraction_VIIRS_NIR_OFOVa" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "cloud_fraction_VIIRS_NIR_OFOVb" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "cloud_fraction_VIIRS_NIR_OFOVc" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "cloud_fraction_VIIRS_SWIR_IFOV" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "cloud_fraction_VIIRS_SWIR_OFOVa" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "cloud_fraction_VIIRS_SWIR_OFOVb" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "cloud_fraction_VIIRS_SWIR_OFOVc" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "eastward_wind" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA", # This apprears to exist in OFFL files but not in RPRO files
        "dry_air_subcolumns" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "fluorescence_apriori" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "instrument_configuration_identifier" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "instrument_configuration_version" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "methane_profile_apriori" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "methane_ratio_weak_strong_standard_deviation" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "methane_strong_twoband_total_column" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "methane_weak_twoband_total_column" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "northward_wind" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA", # This apprears to exist in OFFL files but not in RPRO files
        "pressure_interval" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "reflectance_cirrus_VIIRS_NIR" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "reflectance_cirrus_VIIRS_SWIR" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "scaled_small_pixel_variance" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "surface_altitude" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "surface_altitude_precision" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "surface_classification" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "surface_pressure" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "water_weak_twoband_total_column" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "water_strong_twoband_total_column" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",
        "water_ratio_weak_strong_standard_deviation" : "/PRODUCT/SUPPORT_DATA/INPUT_DATA",

        "methane_mixing_ratio_histogram_axis" : "/METADATA/QA_STATISTICS",
        "methane_mixing_ratio_pdf_axis" : "/METADATA/QA_STATISTICS",
        "methane_mixing_ratio_histogram_bounds" : "/METADATA/QA_STATISTICS",
        "methane_mixing_ratio_pdf_bounds" : "/METADATA/QA_STATISTICS",
        "methane_mixing_ratio_histogram" : "/METADATA/QA_STATISTICS",
        "methane_mixing_ratio_pdf(XCH4_pdf_axis)" : "/METADATA/QA_STATISTICS",

        # /METADATA/ALGORITHM_SETTINGS
        # /METADATA/GRANULE_DESCRIPTION
        # /METADATA/ISO_METADATA
        # /METADATA/EOP_METADATA
        # /METADATA/ESA_METADATA
    }

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.filename = filepath[-86:]
        self.orbit = self.filename[52:57]
        if (not(logger is None)): logger.info("Orbit: {} Processor: {}".format(TROPOMI.get_orbit(self.filename), TROPOMI.get_processor_version(self.filename)))
        ds  = Dataset(filepath)

        self.time_reference_seconds_since_1970 = ds.time_reference_seconds_since_1970
        self.scans = ds[TROPOMI.SP5CH4variables["scanline"]].dimensions["scanline"].size
        self.pixels = ds[TROPOMI.SP5CH4variables["ground_pixel"]].dimensions["ground_pixel"].size
        self.processor_version = TROPOMI.get_processor_version(self.filename)
        ds.close()        
        
        __licencse__ = "GPL"
    
    def chart_box(self, source: Source, geometry: Geometry, boxfilter: Polygon, filename: str) -> None:
        '''
        Chart square box around source. If file

        Parameters
        ----------
        source: Source
            source of emissions
        geometry: Geometry
            Geometry to be used
        size: float
            Size in degrees of the box around the source
        figure: int
            Figure numner
        fileimage: str        
                Destination of image. If None the image is displayed

        '''
        title = "CH4 Mixing Ratio Bias Corrected\nOrbit: {} Processor: {}".format(self.orbit, self.processor_version)

        [minscan, minpixel, maxscan, maxpixel, _, _, lons, lats, scan] = self.narrow_to_domain(geometry, boxfilter)
        
        if (scan is None):
            print("Error: This orbit does not intersects with domain.")
            return
        
        if not (logger is None): logger.info("minscan: {}, maxscan: {}, minpixel:{}, maxpixel: {}".format(minscan, maxscan, minpixel, maxpixel))

        if (0 == np.ma.count(scan)):
            print("Error: All pixels within domain are masked.")
            return

        v_average = np.ma.average(scan)
        if not (logger is None): logger.info(v_average)
        v_min = v_average - 15
        v_max = v_average + 15
        t_lons = boxfilter.exterior.coords.xy[0]
        t_lats = boxfilter.exterior.coords.xy[1]

        e_left = min(t_lons)
        e_right = max(t_lons)
        e_bottom = min(t_lats)
        e_top = max(t_lats)
        lonwidth = e_right - e_left
        latheight = e_top - e_bottom
        mll = lonwidth
        if mll < latheight: mll = latheight
        fracy = int(10 * latheight/mll)
        if fracy < 4: fracy = 4

        plateCarree = ccrs.PlateCarree()
        fig = plt.figure(figsize=(10, fracy))
        ax = plt.axes((0.1, 0, 0.8, 0.9), projection=plateCarree)
        ax.set_extent(extents=[e_left, e_right, e_bottom, e_top], crs=plateCarree)
        ax.coastlines()
        cs = ax.pcolormesh(lons, lats, scan, shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.RdYlBu_r, transform=plateCarree)
        cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
        cbar.ax.set_xlabel("CH4 ppb")

        if boxfilter.contains(source.xy): source.plot_source(ax)
        s_m = create_source("MoranbahNorth")
        if boxfilter.contains(s_m.xy): s_m.plot_source(ax)
        aws = create_aws("Moranbah")
        if boxfilter.contains(aws.xy): aws.plot_aws(ax)

        plt.title(title)
        gl = ax.gridlines(crs=plateCarree, draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')

        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'} 

        if not filename is None:
            plt.savefig(filename)
            print("Chart generated: {}".format(filename))
        else:
            plt.show()
        plt.close("all")
        return

    def chart_domain(self, source: Source, geometry: Geometry, filename:str) -> None:
        '''
        Chart config domain, specifically this reproduces chart Figure 1 subplot f from Sadavarte 2021

        Parameters
        ----------
        source: Source
             
        geometry: Geometry
            Gemoetry used
        figure: int
            Figure numner
        filename: str        
                Destination of image. If None the image is displayed

        '''
        config = Config()
        
        title = "CH4 Mixing Ratio Bias Corrected\nOrbit: {} Processor: {}".format(self.orbit, self.processor_version)
        extent = {"t_lat" : config.Algorithm_CSF_domain_ymax, "l_lon" : config.Algorithm_CSF_domain_xmin, "b_lat" : config.Algorithm_CSF_domain_ymin, "r_lon" : config.Algorithm_CSF_domain_xmax}

        [_, _, _, _, _, _, lons, lats, scan] = self.narrow_to_domain(geometry, config.Algorithm_CSF_domain) # box(146.0, -24.0, 150.0, -20.0))
        if (scan is None):
            print("Error: This orbit does not intersects with domain.")
            return
        
        lonwidth = config.Algorithm_CSF_domain_xmax - config.Algorithm_CSF_domain_xmin
        latheight = config.Algorithm_CSF_domain_ymax - config.Algorithm_CSF_domain_ymin
        mll = lonwidth
        if mll < latheight: mll = latheight
        fracy = int(10 * latheight/mll)
        if fracy < 4: fracy = 4

        plateCarree = ccrs.PlateCarree()
        fig = plt.figure(figsize=(10, fracy))
        ax = plt.axes((0.1, 0, 0.8, 0.9), projection=plateCarree)
        ax.set_extent(extents=[extent["l_lon"], extent["r_lon"], extent["b_lat"], extent["t_lat"]], crs=plateCarree)
        ax.coastlines()
        cs = ax.pcolormesh(lons, lats, scan, shading="nearest", vmin=1793, vmax=1828, cmap=plt.cm.RdYlBu_r, transform=plateCarree)
        cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
        cbar.ax.set_xlabel("CH4 ppb")

        if config.Algorithm_CSF_domain.contains(source.xy): source.plot_source(ax)
        s_m = create_source("MoranbahNorth")
        if config.Algorithm_CSF_domain.contains(s_m.xy): s_m.plot_source(ax)
        aws = create_aws("Moranbah")
        if config.Algorithm_CSF_domain.contains(aws.xy): aws.plot_aws(ax)

        plt.title(title)
        gl = ax.gridlines(crs=plateCarree, draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')

        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'} 

        if not filename is None:
            plt.savefig(filename)
            print("Chart generated: {}".format(filename))
        else:
            plt.show()
        plt.close("all") 
        return 

    def chart_size(self, source: Source, geometry: Geometry, size: float, filename: str) -> None:
        '''
        Chart square box around source using size in degrees

        Parameters
        ----------
        source: Source
            source of emissions
        geometry: Geometry
            Geometry to be used
        size: float
            Size in degrees of the box around the source
        figure: int
            Figure numner
        filename: str        
                Destination of image. If None the image is displayed

        '''
        title = "CH4 Mixing Ratio Bias Corrected\nOrbit: {} Processor: {}".format(self.orbit, self.processor_version)

        [minscan, minpixel, maxscan, maxpixel, _, _, lons, lats, scan] = self.narrow_to_domain(geometry, box(source.xy.x - size, source.xy.y - size,  source.xy.x + size, source.xy.y + size))
        
        if (scan is None):
            print("Error: This orbit does not intersects with domain.")
            return
        
        if not (logger is None): logger.info("minscan: {}, maxscan: {}, minpixel:{}, maxpixel: {}".format(minscan, maxscan, minpixel, maxpixel))

        if (0 == np.ma.count(scan)):
            print("Error: All pixels within domain are masked.")
            return

        v_average = np.ma.average(scan)
        v_min = v_average - 15
        v_max = v_average + 15

        plateCarree = ccrs.PlateCarree()
        fig = plt.figure(figsize=(10,8))
        ax = plt.axes((0.1, 0, 0.8, 0.9), projection=plateCarree,)
        ax.set_extent(extents=[source.xy.x - size, source.xy.x + size , source.xy.y - size, source.xy.y + size], crs=plateCarree)
        ax.coastlines()
        cs = ax.pcolormesh(lons, lats, scan, shading="nearest", vmin=v_min, vmax=v_max, cmap=plt.cm.RdYlBu_r, transform=plateCarree)
        cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
        cbar.ax.set_xlabel("CH4 ppb")

        poly = box(source.xy.x - size, source.xy.y - size, source.xy.x + size, source.xy.y + size)
        if poly.contains(source.xy): source.plot_source(ax)

        s_m = create_source("MoranbahNorth")
        if poly.contains(s_m.xy): s_m.plot_source(ax)

        aws = create_aws("Moranbah")
        if poly.contains(aws.xy): aws.plot_aws(ax)

        plt.title(title)
        gl = ax.gridlines(crs=plateCarree, draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')

        gl.top_labels = False
        gl.left_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10, 'color': 'gray'} 

        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)
            print("Chart generated: {}".format(filename))
        plt.close("all")        
        return 

    def describe_path(self, p: str) -> None:
        '''
        Print description of netCDF4-d variable

        Parameters
        ----------
        p:str
            Path of the element

        Returns
        -------
        None

        '''
        try:
            ds  = Dataset(self.filepath) 
            if (p): print(ds[p])
            else: print(ds)
        except IndexError:
            print("Error: Path {} does not exist.".format(p))
        finally:
            if not(ds is None): ds.close()
        pass

    def describe_variable(self, variablename) -> None:
        '''
        Print description of netCDF4-d variable

        Parameters
        ----------
        variablename:str
            One of the elements of TROPOMI.SP5CH4variables

        Returns
        -------
        None

        '''
        try:
            ds  = Dataset(self.filepath) 
            print(ds[TROPOMI.SP5CH4variables[variablename]].variables[variablename])
        finally:
            if not(ds is None): ds.close()

    def get_pixel_for_source(self, source: Point) -> tuple[int, int, datetime, datetime]:
        '''
        Find pixel containing source
        
        Parameters
        ----------
        source: Point
            Point for which we are coordinates we are looking for

        Returns
        -------
        tuple[scan, pixel, datetime, datetime]
            - scan: int Scan coordinate
            - pixel: int Pixel coordinate 
            - time_at_source: datetime Time of pixel at source
            - time_at_source_nearest_hour: datetime Time of pixel at source rounded to the nearest hour
        '''
        s = 0
        lt_box_arr = self.get_variable_data("latitude_bounds")[0, :, :, :]
        ln_box_arr = self.get_variable_data("longitude_bounds")[0, :, :, :]
        delta_time = self.get_variable_data("delta_time")[0,:]
        mxlt = np.max(lt_box_arr)
        mnlt = np.min(lt_box_arr)
        mxln = np.max(ln_box_arr)
        mnlon = np.min(ln_box_arr)
        
        if mxlt < source.y: return None
        if source.y < mnlt: return None
        if mxln < source.x: return None
        if source.x < mnlon: return None
        
        while s < self.scans:
            p = 0
            while p < self.pixels:
                lats = lt_box_arr[s,p,:]
                lons = ln_box_arr[s,p,:]
                mxlt = max(lats) 
                mnlt = min(lats)
                mxln = max(lons)
                mnlon = min(lons)
                
                if mxlt < source.y or source.y < mnlt or mxln < source.x or  source.x < mnlon: 
                    p += 1
                    continue

                # See file:///Z:/Sentinel5/Documentation/Sentinel-5P-Level-2-Input-Output-Data-Definition.pdf page 141
                bl = (lons[0], lats[0])
                br = (lons[1], lats[1])
                ur = (lons[2], lats[2])
                ul = (lons[3], lats[3])
                pixelbox = Polygon([bl,br,ur,ul])
                if  (pixelbox.contains(source)):
                    delta = int(delta_time[s])
                    dt_exact = datetime.fromtimestamp(self.time_reference_seconds_since_1970, tz=timezone.utc) + timedelta(milliseconds=delta)
                    dt_rounded = datetime(dt_exact.year, dt_exact.month, dt_exact.day, dt_exact.hour, tzinfo=timezone.utc)
                    if 30 <= dt_exact.minute: dt_rounded += timedelta(hours=1)
                    return s, p, dt_exact, dt_rounded
                p += 1
            s += 1
        return None, None, None, None

    def get_variable_data(self, variablename) -> np.array:
        '''
        Get netCDF4-d variable

        Parameters
        ----------
        variablename:str
            One of the elements of ROPOMI.SP5CH4variables

        Returns
        -------
        data: numpy. array
            netCDF4 variable data

        '''
        ds = None
        try:
            data = None
            ds  = Dataset(self.filepath) 
            data = ds[TROPOMI.SP5CH4variables[variablename]].variables[variablename][:]
        finally:
            if not(ds is None): ds.close()
        return data
        
    def narrow_to_domain(self, geometry: Geometry, domainbox: MultiPolygon) -> tuple[int, int, int, int, np.array, np.array, np.array, np.array, np.ma.array]: 
        '''
        Get pixels in a box defined by domainbox for TROPOMI variable from a file 

        Parameters
        ----------
        geometry (Geometry)
            - center - inlucde pixels with centre within domain only
            - contains - inlucde pixels conained within domain only
            - intersects - include pixels intersecting with domain
        domainbox (MultiPolygon) - MultiPolygon limiting pixels

        Returns
        -------
        tuple containing
            - minscan int: minimum scan dimension inclusive
            - minpixel int: minimum pixel dimension inclusive
            - maxscan int: maximum scan dimension exlusive
            - maxpixel int: maximum pixel dimesnion exlusive
            - lonbox (array[minscan:maxscan, minpixel:maxpixel, 4]): Longitudes of pixels corners
            - latbox (array[minscan:maxscan, minpixel:maxpixel, 4]): Latitudes of pixels corners
            - loncen (array[minscan:maxscan, minpixel:maxpixel]): Longitudes of pixels centers
            - latcen(array[minscan:maxscan, minpixel:maxpixel]): Latitudes of pixels centers
            - scan (numpy.ma.array[minscan:maxscan, minpixel:maxpixel]) - Masked array of variable pixels limited to domainbox using geometry
        '''
        prepare(domainbox)
        minscan = self.scans + 1
        maxscan = 0
        minpixel = self.pixels + 1
        maxpixel = 0

        variable = self.get_variable_data("methane_mixing_ratio_bias_corrected")[0, :, :]
        lt_box_arr = self.get_variable_data("latitude_bounds")[0, :, :, :]
        ln_box_arr = self.get_variable_data("longitude_bounds")[0, :, :, :]
        ln_center_arr = self.get_variable_data("longitude")[0, :, :]
        lt_center_arr = self.get_variable_data("latitude")[0, :, :]              

        [minlon, minlat, maxlon, maxlat] = domainbox.bounds
        for s in range(self.scans):
            
            # Optimize
            if np.max(lt_box_arr[s,:,:]) < minlat: continue
            if maxlat < np.min(lt_box_arr[s,:,:]): continue
            
            for p in range(self.pixels):
                # Optimize
                if np.max(ln_box_arr[s,p,:]) < minlon: continue
                if maxlon <  np.min(ln_box_arr[s,p,:]): continue

                # See file:///Z:/Sentinel5/Documentation/Sentinel-5P-Level-2-Input-Output-Data-Definition.pdf page 141
                ul = (ln_box_arr[s,p,3], lt_box_arr[s,p,3])
                bl = (ln_box_arr[s,p,0], lt_box_arr[s,p,0])
                br = (ln_box_arr[s,p,1], lt_box_arr[s,p,1])
                ur = (ln_box_arr[s,p,2], lt_box_arr[s,p,2])
                pixelbox = Polygon([ul,bl,br,ur])
                center = Point(ln_center_arr[s,p], lt_center_arr[s,p])
                if ((geometry == Geometry.CENTER) and not(domainbox.contains(center))): 
                    variable[s,p] = np.ma.masked
                    continue
                if ((geometry == Geometry.INTERSECTS) and not(domainbox.intersects(pixelbox))): 
                    variable[s,p] = np.ma.masked
                    continue
                if ((geometry == Geometry.CONTAINS) and not(domainbox.contains(pixelbox))): 
                    variable[s,p] = np.ma.masked
                    continue
                if p < minpixel : minpixel = p
                if p >= maxpixel: maxpixel = p + 1
                if s < minscan: minscan = s
                if s >= maxscan: maxscan = s + 1

        if (0 == maxscan or 0 == maxpixel): 
            if (not(logger is None)): 
                logger.error("No data for geometry: {}".format(geometry.name))
                if domainbox is MultiPolygon:
                    for poly in domainbox:
                        points = poly.exterior.coords[:-1]
                        logger.error("Polygon")
                        for p in points:
                            logger.error("\tPoint: [{}, {}]".format(p.x, p.y))
                elif domainbox is Polygon:
                    points = domainbox.exterior.coords[:-1]
                    logger.error("Polygon")
                    for p in points:
                        logger.error("\tPoint: [{}, {}]".format(p.x, p.y))
            return None, None, None, None, None, None, None, None, None

        values = variable[minscan : maxscan, minpixel : maxpixel]
        latcen = lt_center_arr[minscan : maxscan, minpixel : maxpixel]
        loncen = ln_center_arr[minscan : maxscan, minpixel : maxpixel]
        latbox = lt_box_arr[minscan : maxscan, minpixel : maxpixel, :]
        lonbox = ln_box_arr[minscan : maxscan, minpixel : maxpixel, :]

        maxvalues = np.ma.max(values)
        minvalues = np.ma.min(values)

        if maxvalues == np.ma.masked: smax = "--"
        elif np.isfinite(maxvalues): smax = "{:.3f}".format(maxvalues) 
        else: smax = "{}".format(maxvalues)

        if minvalues == np.ma.masked: smin = "--"
        elif np.isfinite(minvalues): smin = "{:.3f}".format(minvalues) 
        else: smin = "{}".format(minvalues)

        message = "{} stats: scan={}:{}, pixel={}:{}, range={}:{}".format(geometry, minscan, maxscan, minpixel, maxpixel, smin, smax)
        if (not(logger is None)): logger.info(message)

        return minscan, minpixel, maxscan, maxpixel, lonbox, latbox, loncen, latcen, values

    @staticmethod
    def exists_locally(config: Config, orbit:str, processorversion:str) -> bool:
        '''
        Check if this file exists locally
        
        Parameters
        ----------
        orbit: str
            Orbit as a 5 digit string
        processorversion: str
            If None any version will do. If specfied use specific version in format e.g. 010301 MMmmrr - M major version, m minor version, r release 
        '''
        folder = path.join(Config.TROPOMI_folder, processorversion)
        if not path.exists(folder): return False
        filenames = listdir(folder)
        for f in filenames:
            if not len(f) == 86: continue
            if not TROPOMI.get_extension(f) == "nc": continue
            if not TROPOMI.get_productidentifier(f) == "L2__CH4____": continue
            if TROPOMI.get_orbit(f) == orbit:
                if processorversion is None: return True
                elif TROPOMI.get_processor_version(f) == processorversion: return True
        
        return False
    
    @staticmethod
    def get_collection_number(filename: str) -> str:
        '''
        Get collection number version time from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''        
        return filename[58:60]
    
    @staticmethod
    def get_extension(filename: str) -> str:
        '''
        Get extension from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''        
        return filename[84:86]

    @staticmethod
    def get_missionname(filename: str) -> str:
        '''
        Get mission name version time from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''        
        return filename[0:3]
    
    @staticmethod
    def get_orbit(filename: str) -> str:
        '''        
        Get orbit from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''
        return filename[52:57]

    @staticmethod
    def get_orbit_int(filename: str) -> int:
        '''        
        Get orbit from file name as integer
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''
        return int(filename[52:57])

    @staticmethod
    def get_processing_stream(filename: str) -> str:
        '''
        Get processing stream NRTI (near real time) OFFL (offline)  RPRO (reprocessing) from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''        
        return filename[4:8]
    
    @staticmethod
    def get_processing_time(filename: str) -> str:
        '''
        Get processing time from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''        
        return filename[68:83]

    @staticmethod
    def get_processor_version(filename: str) -> str:
        '''
        Get processor version time from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''        
        return filename[61:67]

    @staticmethod
    def get_productidentifier(filename: str) -> str:
        '''        
        Get product identifier from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11

        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''
        return filename[9:20]

    @staticmethod
    def get_scan_end(filename: str) -> str:
        '''        
        Get orbit end time from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''
        return filename[36:51]

    @staticmethod
    def get_scan_start(filename: str) -> str:
        '''        
        Get orbit start time from file name
        Sentinel-5P-Level-2-Product-User-Manual-Methane.pdf page 11
        0         1         2         3         4         5         6         7         8
        01234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        S5P_OFFL_L2__CH4____20190915T030716_20190915T044845_09956_01_010302_20190921T050748.nc
        '''
        return filename[20:35]

def TROPOMI_for_orbit(orbit:str, processor: str) -> TROPOMI:
    '''
    Factory method returning TROPOMI using orbit.

    Parameters
    ----------
    orbit: str
        Orbit number as a 5 digit string. e.g. 09956
    processor: str
        Processor number as a 6 digit string e.g. 020400

    Returns
    -------
    TOPOMI class or raises an exception
    '''

    # Try to find the file in directory
    filepath = None
    sourcedir = path.join(Config.TROPOMI_folder, processor)
    for f in listdir(sourcedir):
        if len(f) == 86 and f.startswith("S5P") and f.endswith(".nc") and 0 < f.index("L2__CH4___"):
            if TROPOMI.get_orbit(f) == orbit:
                filepath = path.join(sourcedir, f)
                if not(logger is None): logger.info("File found: {}".format(f))
                break

    if filepath is None:
        if not(logger is None): logger.error("No TROPOMI file found orbit: {}".format(orbit))
        raise RuntimeError("No TROPOMI file found orbit: {}".format(orbit))

    return TROPOMI(filepath)

def _handler_chart_box(args: Namespace):
    '''
    Handler for TROPOMI charting

    Parameters
    ----------
    args (dict) - List of arguments
    '''
    source = create_source(args.source)
    if args.geometry is None: geometry = Geometry.INTERSECTS
    if args.geometry == "center": geometry = Geometry.CENTER
    if args.geometry == "contains": geometry = Geometry.CONTAINS
    if args.geometry == "intersects": geometry = Geometry.INTERSECTS
    
    Config.create_log( "TROPOMI_" + source.case_name +"_" + args.orbit +".log")
    boxfilter = box(args.lb_lon, args.lb_lat, args.ut_lon, args.ut_lat)
    tropomi = TROPOMI_for_orbit(args.orbit, args.processor)
    tropomi.chart_box(source, geometry, boxfilter, None)

def _handler_chart_domain(args: Namespace):
    '''
    Handler for TROPOMI charting

    Parameters
    ----------
    args (dict) - List of arguments
    '''
    source = create_source(args.source)
    if args.geometry is None: geometry = Geometry.INTERSECTS
    if args.geometry == "center": geometry = Geometry.CENTER
    if args.geometry == "contains": geometry = Geometry.CONTAINS
    if args.geometry == "intersects": geometry = Geometry.INTERSECTS

    Config.create_log("TROPOMI_" + source.case_name +"_" + args.orbit +".log")
    tropomi = TROPOMI_for_orbit(args.orbit, args.processor)
    tropomi.chart_domain(source, geometry, None)

def _handler_chart_size(args: Namespace):
    '''
    Handler for TROPOMI charting

    Parameters
    ----------
    args (dict) - List of arguments
    '''
    source = create_source(args.source)
    if args.geometry is None: geometry = Geometry.INTERSECTS
    if args.geometry == "center": geometry = Geometry.CENTER
    if args.geometry == "contains": geometry = Geometry.CONTAINS
    if args.geometry == "intersects": geometry = Geometry.INTERSECTS

    Config.create_log("TROPOMI_" + source.case_name +"_" + args.orbit +".log")
    tropomi = TROPOMI_for_orbit(args.orbit, args.processor)
    tropomi.chart_size(source, geometry, args.size, None)

if __name__ == '__main__':
    parser = ArgumentParser(prog="TROPOMI", description="Chart TROPOMI data")

    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_chart = subparsers.add_parser("chart", help="Chart commands.")
    subparsers_chart = parser_chart.add_subparsers(help="subcommand help", required=True)

    parser_box = subparsers_chart.add_parser("box")
    parser_box.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_box.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 04579")
    parser_box.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_box.add_argument("lb_lon", type=float, help="Lower bottom longitude" )
    parser_box.add_argument("lb_lat", type=float, help="lower bottom latitude" )
    parser_box.add_argument("ut_lon", type=float, help="Upper top longitude" )
    parser_box.add_argument("ut_lat", type=float, help="Upper top latitude" )
    parser_box.add_argument("-g","--geometry", type=str, choices=["center", "contains", "intersects"], help="Geometry to be used for chart")
    parser_box.set_defaults(func=_handler_chart_box)

    parser_domain = subparsers_chart.add_parser("domain")
    parser_domain.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_domain.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 04579")
    parser_domain.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_domain.add_argument("-g","--geometry", type=str, choices=["center", "contains", "intersects"], help="Geometry to be used for chart")
    parser_domain.set_defaults(func=_handler_chart_domain)

    parser_size = subparsers_chart.add_parser("size")
    parser_size.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_size.add_argument("orbit", type=str, help="Orbit number as a five digit string e.g. 04579")
    parser_size.add_argument("processor", type=str, help="Processor number as a six digit string e.g. 020400")
    parser_size.add_argument("size", type=float, help="Size of the box diagram 1 - one degree in each direction, 0.5 - half degree in each direction" )
    parser_size.add_argument("-g","--geometry", type=str, choices=["center", "contains", "intersects"], help="Geometry to be used for chart")
    parser_size.set_defaults(func=_handler_chart_size)

    args = parser.parse_args()

    if not(len(args.orbit) == 5) or not(args.orbit.isdigit()):
        print ("Orbit must be a five digit number: {}".format(args.orbit))
        exit()

    if not(len(args.processor) == 6) or not(args.processor.isdigit()):
        print ("Processor must be a 6 digit number: {}".format(args.processor))
        exit()

    args.func(args)
