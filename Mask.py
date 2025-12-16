
from enum import Enum
class Mask(Enum):
    '''
    Calculation of enhancement after subtraction of background
    1 - Do not mask
    2 - Mask negative values
    '''
    NONE = 1
    NEGATIVE = 2