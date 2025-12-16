from Constants import GRAVITY
import logging
from math import atan2, pi
logger = logging.getLogger(__name__)

class Meteorology:
    @staticmethod
    def __calculate_PressureWeightedAverage(p_bottom: float, p_top: float, p : list[float] , f: list[float]) -> float:
        '''
        Calculate pressure weighted averaged between two pressure levels

        Parameters
        ----------
        p_bottom: float
            Pressure at the lowest included level
        p_top: float
            Pressure at the highest included level
        p: array
            1D - Array of pressure
        f: array
            1D - Array of variable values on pressure levels

        Returns
        ------
        av - average or 0
        '''
        plen = len(p)
        pint = 0
        fint = 0
        for i in range(plen):
            if p_top <= p[i] and p[i] <= p_bottom:
                deltap = 0 if i == 0 or i == (plen - 1) else 0.5 * (p[i - 1] - p[i + 1])
                fint += deltap * f[i] * p[i]
                pint += deltap * p[i]

        return fint / pint if not(pint == 0) else 0
    
    @staticmethod
    def calculate_ERA5_PressureAveragedWind(p: list[float], z: list[float], u: list[float], v: list[float], blh: float, sp:float) -> tuple[float, float]:
        '''
        Calculate pressure averaged wind 

        Parameters
        ----------
        p: list[float]
            List of pressure levels in descending order in hPa
        z: list[float]
            List of geopotential for each pressure level
        u: list[float]
            List of west-east wind component for each pressure level m/s
        v: list[float]
            List of south-north wind component for each pressure level m/s
        blh: float
            Boundary layer height in m
        sp: float
            Surface pressure in hPa

        Returns
        -------
        [u,v]: array
            U and V component of wind
        '''
        if not(logger is None): logger.info("")
        
        era5_p_len = len(p)

        pindex_blh = 0
        
        for i  in range(era5_p_len):
            if z[i] / GRAVITY < blh:
                pindex_blh = i

        ptop = p[pindex_blh]
        
        # Find the pressure of the lowest level above station level pressure
    
        pindex_sp = 0
        while pindex_sp < era5_p_len:
            if sp < p[pindex_sp]:
                pindex_sp += 1
            else:
                break

        pbottom = p[pindex_sp]

        if not(logger is None): logger.info("Pressure bottom: pmax: {} surface_pressure: {} (hPa)".format(pbottom, sp))
        if not(logger is None): logger.info("Pressure top:    pmax: {} blh:{} (m)".format(ptop, blh))

        for i in range(pindex_sp, pindex_blh + 1):
            if not(logger is None): logger.info("Pressure: {} (hPa) Height: {} (m): U: {} (m/s) V: {} (m/s)".format(p[i], z[i] / GRAVITY, u[i], v[i]))

        uav = Meteorology.__calculate_PressureWeightedAverage(pbottom, ptop, p, u)
        vav = Meteorology.__calculate_PressureWeightedAverage(pbottom, ptop, p, v)
        if not(logger is None): logger.info("Pressure Weighted Average U: {} (m/s) V: {} (m/s)".format(uav, vav))
        return [uav, vav]

    def calculate_azimuth_degree_meteorology(u: float, v: float) -> float:
        '''
        Calculate azimuth in degrees using meteorological convention
        0 - represents wind blowing from the north u = 0, v = -1
        90 - represents wind blowing from the east u = -1, v = 0
        180 - represents wind blowing form the south u = 0, v = 1
        270 - represents wind blowing form the west u = 1, v = 0

        Parameters
        ----------
        u: float
            West to East wind component
        v: float
            South to North wind component

        Returns
        -------
        azimuth: 
        '''
        az = atan2(u,v) * 180. / pi + 180.0 # Meteorological convention
        if az < 0: az += 360.0
        if 360.0 <= az: az -= 360.0
        return az
    
    def calculate_azimuth_degree_oceanography(u: float, v: float) -> float:
        '''
        Calculate azimuth in degrees using oceanographic convention
        0 - represents wind blowing to the north u = 0, v = 1
        90 - represents wind blowing to the east u = 1, v = 0
        180 - represents wind blowing to the south u = 0, v = -1
        270 - represents wind blowing to the west u = -1, v = 0

        Parameters
        ----------
        u: float
            West to East wind component
        v: float
            South to North wind component

        Returns
        -------
        azimuth: 
        '''
        az = atan2(u,v) * 180. / pi  # Oceanogrphic convention
        if az < 0: az += 360.0
        if 360.0 <= az: az -= 360.0
        return az