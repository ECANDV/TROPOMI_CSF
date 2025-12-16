from Constants import CH4_MOLECULAR_MASS
from Constants import DRYAIR_MOLECULAR_MASS
from Constants import H2O_MOLECULAR_MASS
from Constants import GRAVITY

def convert_column_ppb_dry_air(ppb: float, p: float) -> float:
    '''
    Convert column ppb enhancement to concentration enhancement in kg m-2

    Parameters
    ----------
    ppb: float
        Concentration in ppb
    p: float
        Pressure in Pa (N m-2 = kg m s-2 m-2 = kg m-1 s-2 )

    Returns
    -------
    Mass of methane in unit area column (kg m-2) = (kg m-1 s-2) / (m s-2)
    '''
    return ppb * CH4_MOLECULAR_MASS / DRYAIR_MOLECULAR_MASS * p / GRAVITY * 1E-9

def convert_column_ppb_with_water(ppb: float, dry_air: float, h2o:float, p: float) -> float:
    '''
    Convert column ppb enhancement to concentration enhancement in kg m-2

    Parameters
    ----------
    ppb: float
        Concentration in ppb
    dry_air: float
        Dry air column in (mol m^-2) as average of TROPOMI "dry_air_subcolumns"
    h2O: float
        Water column in (mol m^-2) from TROPOMI "water_total_column"
    p: float
        Pressure in Pa (N m-2 = kg m s-2 m-2 = kg m-1 s-2 )

    Returns
    -------
    Mass of methane in unit area column (kg m-2) = (kg m-1 s-2) / (m s-2)
    '''
    return ppb * CH4_MOLECULAR_MASS * ((dry_air + h2o) / (dry_air * DRYAIR_MOLECULAR_MASS + h2o * H2O_MOLECULAR_MASS)) * p / GRAVITY * 1E-9
