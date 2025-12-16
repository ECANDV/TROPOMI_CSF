from argparse import ArgumentParser, Namespace
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import logging
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import numpy as np
from os import path
from osgeo import gdal

logger = logging.getLogger(__name__)

path_9s_geotiff = path.join("F:","Auslig", "DEM-9S", "Data_9secDEM_D8", "dem-9s","dem-9s.tif")

class SRTM:
    '''
    Shuttle Radar Topgraphy Mission DEM data download from auslig
    '''
    def __init__(self, t_lat, l_lon, b_lat, r_lon, elevation):
        self.elevation = elevation
        self.t_lat = t_lat
        self.l_lon = l_lon
        self.b_lat = b_lat
        self.r_lon = r_lon

    def chart(self, fig: plt.Figure, ax: Axes):
        plateCarree = ccrs.PlateCarree()
        ax.set_extent(extents=(self.l_lon, self.r_lon, self.b_lat, self.t_lat), crs=plateCarree)
        # use where() to extract 2-dimensional subsets of elev_arr
        ele = np.where(0 < self.elevation, self.elevation, np.nan)
        if not logger is None: logger.info("\tMinimum: {} (m) Maximum: {} (m)".format(np.min(ele), np.max(ele)))

        cs = ax.imshow(ele, cmap=plt.cm.Blues, transform=plateCarree, extent=(self.l_lon, self.r_lon, self.b_lat, self.t_lat))
        cbar = fig.colorbar(cs, orientation="horizontal", pad=0.075)
        cbar.ax.set_xlabel("m")

def SRTM_from_geotiff(filename, t_lat, l_lon, b_lat, r_lon) -> SRTM:
    '''
    Factory method for SRTM
    '''
    if not logger is None: logger.info("GeoTIFF: {}".format(filename))
    if not path.exists(filename):
        print("Error: missing file {}".format(filename))
        if not logger is None: logger.error("Error: missing file {}".format(filename))
        return None

    gdal.UseExceptions()
    # Get dataset
    ds = gdal.Open(filename,  gdal.GA_ReadOnly)
    if (None == ds):
        print("\tDataset is empty")
        return None

    size_lon = ds.RasterXSize
    size_lat = ds.RasterYSize
    if not logger is None: logger.info("\tTotal size is {} x {} x {}".format(size_lon, size_lat, ds.RasterCount))

    geotransform = ds.GetGeoTransform()
    ul_lon = geotransform[0]
    ul_lat = geotransform[3]
    px_lon = geotransform[1]
    px_lat = geotransform[5]
    pad_lon = geotransform[2]
    pad_lat = geotransform[4]

    px_x_min = np.iinfo(np.int32).max
    px_x_max = -1
    px_y_min = np.iinfo(np.int32).max
    px_y_max = -1

    if geotransform:
        if not logger is None: logger.info("\tOrigin = ({}, {})".format(ul_lon, ul_lat))
        if not logger is None: logger.info("\tPixel Size = ({}, {})".format(px_lon, px_lat))
        if not logger is None: logger.info("\tTransformation = ({}, {})".format(pad_lon, pad_lat))

    band = ds.GetRasterBand(1)

    # In this case pad_lon and pad_lat are 0 so we can do this independently
    for n in range(size_lon):
        lon = ul_lon + n * px_lon
        if l_lon <= lon and lon <= r_lon:
            if (n < px_x_min):
                px_x_min = n
            if (px_x_max < n):
                px_x_max = n

    for t in range (size_lat):
        lat = ul_lat + t * px_lat
        if b_lat <= lat and lat <= t_lat:
            if (t < px_y_min):
                px_y_min = t
            if (px_y_max < t):
                px_y_max = t

    w_xsize = px_x_max - px_x_min + 1
    w_ysize = px_y_max - px_y_min + 1

    if not logger is None: logger.info("\tBox x_min: {} x_max: {} y_min: {} y_max: {}".format(px_x_min, px_x_max, px_y_min, px_y_max))
    elevation = band.ReadAsArray(
        xoff = px_x_min, 
        yoff = px_y_min, 
        win_xsize = w_xsize, 
        win_ysize = w_ysize,
        buf_xsize = w_xsize,
        buf_ysize = w_ysize,
        buf_type = gdal.GDT_Float64
    )
    ds.Close()
    return SRTM(t_lat, l_lon, b_lat, r_lon, elevation)

def _handler_chart(args: Namespace):
    srtm = SRTM_from_geotiff(path_9s_geotiff, -20.6, 146, -24, 150)
    fig = plt.figure(figsize=(10, 10))
    plateCarree =  ccrs.PlateCarree()
    ax = plt.axes(projection=plateCarree)
    srtm.chart(fig, ax)
    plt.show()
    plt.close("all")
    pass

if __name__ == '__main__':
    '''
    This is just for smilarity with other objects
    '''
    parser = ArgumentParser(prog="SRTM")
    parser.set_defaults(func=_handler_chart)
    args = parser.parse_args()
    args.func(args)
