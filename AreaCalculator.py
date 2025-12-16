import cartopy.crs as ccrs
import math
import pyproj
from shapely.geometry import Polygon
from Source import Source

class AreaCalculator:
    '''
    Calculate area in m^2 of of a polygon expressed in lon/lat geodetic coordinates.
    UTM is the area preserving coordinate so convert to UTM and then use shapely to get the area.
    '''
    @staticmethod
    def calculate_area(source: Source, original_polygon: Polygon ) -> float:

        utm = math.ceil((source.longitude + 180) / 6)
        target_crs = ccrs.UTM(zone=utm)
        original_crs = ccrs.Geodetic()

        # 3. Transform the polygon
        project = pyproj.Transformer.from_crs(original_crs, target_crs, always_xy=True)
        (px, py) = project.transform(original_polygon.exterior.coords.xy[0] , original_polygon.exterior.coords.xy[1])
        projected_polygon = Polygon(list(zip(px, py)))

        return projected_polygon.area