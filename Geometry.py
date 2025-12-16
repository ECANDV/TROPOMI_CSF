from enum import Enum
class Geometry(Enum):
    '''
    Geometrical relationship between domain and pixels
    1 - Pixel as a polygon and is fully contained within domain
    2 - Center of a pixel is within domain
    3 - Pixel as a polygon intersects with domain
    '''
    CONTAINS = 1
    CENTER = 2
    INTERSECTS = 3
