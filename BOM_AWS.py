from argparse import ArgumentParser, Namespace
from Config import Config
import csv
from datetime import datetime, timedelta, timezone
import logging
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import inf, linspace, pi, zeros
from shapely import Point
from os import path

logger = logging.getLogger(__name__)

# Moranbah AWS data can be purchased from Australian Bureau of Meteorology
path_Moranbah_AWS_CSV = path.join("F:", "BoM", "MORANBAH_AIRPORT_20180101_20210729.csv")

# Bins for rose wind
spd_bins = [0, 10, 20, 30, inf]

class AWS_Observation:
    '''
    Single observation based on Moranbah AWS data
    '''
    def __init__(self, displayname: str, date: datetime, longitude: float, latitude: float):
        '''
        Constructor
        '''
        self.displayname = displayname
        self.longitude = longitude
        self.latitude = latitude
        self.Ceilometer_1_amount: str = None
        self.Ceilometer_1_height: int = None # m
        self.Ceilometer_2_amount: str = None
        self.Ceilometer_2_height: int = None # m
        self.Ceilometer_3_amount: str = None
        self.Ceilometer_3_height: int = None # m
        self.Precip_10minutes: float = None # mm
        self.Precip_since9am: float = None # mm
        self.Pressure_MSLP: str = None # hPa
        self.Pressure_QNH: str = None # hPa
        self.Pressure_Station: str = None # hPa
        self.RelHum: int = None # %
        self.Temperature_Air: float = None # C
        self.Temperature_DewPoint: float = None # C
        self.Temperature_WetBulb: float = None # c
        self.UTC_Datetime: datetime = date
        self.Vis: float = None # km
        self.Wind_Gust: float = None # km/h
        self.Wind_Direction_Deg: int = None
        self.Wind_Speed: float = None # km/h

    def __str__(self):
        m = []
        if not self.Temperature_Air is None: m.append("\tTemperature Air: {} (C)".format(self.Temperature_Air))
        if not self.Temperature_WetBulb is None: m.append("\tTemperature Wet Bulb: {} (C)".format(self.Temperature_WetBulb))
        if not self.Temperature_DewPoint is None: m.append("\tDew Point: {} (C)".format(self.Temperature_DewPoint))
        if not self.RelHum is None: m.append("\tRelative Humidity: {} (%)".format(self.RelHum))

        if not self.Precip_10minutes is None: m.append("\tPrecipitation 10min: {} (mm)".format(self.Precip_10minutes))
        if not self.Precip_since9am is None: m.append("\tPrecipitation since 9am (Local): {} (mm)".format(self.Precip_since9am))

        if not self.Ceilometer_1_amount: m.append("\tCloud amount (group 1): {}".format(self.Ceilometer_1_amount))
        if not self.Ceilometer_1_height is None: m.append("\tCloud height (group 1): {} (m)".format(self.Ceilometer_1_height))

        if not self.Ceilometer_2_amount: m.append("\tCloud amount (group 2): {}".format(self.Ceilometer_2_amount))
        if not self.Ceilometer_2_height is None: m.append("\tCloud height (group 2): {} (m)".format(self.Ceilometer_2_height))

        if not self.Ceilometer_3_amount: m.append("\tCloud amount (group 3): {}".format(self.Ceilometer_3_amount))
        if not self.Ceilometer_3_height is None: m.append("\tCloud height (group 3): {} (m)".format(self.Ceilometer_3_height))

        if not self.Wind_Speed is None: m.append("\tWind speed: {} (km/h)".format(self.Wind_Speed))
        if not self.Wind_Direction_Deg is None: m.append("\tWind direction: {} (True deg)".format(self.Wind_Direction_Deg))
        if not self.Wind_Gust is None: m.append("\tMax wind gust: {} (km/h)".format(self.Wind_Gust))

        if not self.Vis is None: m.append("\tVisibility: {} (km)".format(self.Vis))
        else: m.append("\tVisibility: > 10 (km)")

        if not self.Pressure_MSLP is None: m.append("\tPressure Mean Sea Level: {} (hPa)".format(self.Pressure_MSLP))
        if not self.Pressure_Station is None: m.append("\tPressure Station Level: {} (hPa)".format(self.Pressure_Station))
        if not self.Pressure_QNH is None: m.append("\tPressure QNH: {} (hPa)".format(self.Pressure_QNH))
        return "\n".join(m)
    
    def to_csv(self):
        '''
        Default string conversion
        '''
        m = []
        m.append("{} (UTC)".format(self.UTC_Datetime.strftime("%Y-%m-%d %H:%M")))
        
        if not self.Temperature_Air is None: m.append("{} (C)".format(self.Temperature_Air))
        else: m.append("")
        if not self.Temperature_WetBulb is None: m.append("{} (C)".format(self.Temperature_WetBulb))
        else: m.append("")
        if not self.Temperature_DewPoint is None: m.append("{} (C)".format(self.Temperature_DewPoint))
        else: m.append("")
        if not self.RelHum is None: m.append("{} (%)".format(self.RelHum))
        else: m.append("")

        if not self.Precip_10minutes is None: m.append("{} (mm)".format(self.Precip_10minutes))
        else: m.append("")
        if not self.Precip_since9am is None: m.append("{} (mm)".format(self.Precip_since9am))
        else: m.append("")

        if not self.Wind_Speed is None: m.append("{} (km/h)".format(self.Wind_Speed))
        else: m.append("")
        if not self.Wind_Direction_Deg is None: m.append("{} (True deg)".format(self.Wind_Direction_Deg))
        else: m.append("")
        if not self.Wind_Gust is None: m.append("{} (km/h)".format(self.Wind_Gust))
        else: m.append("")

        if not self.Ceilometer_1_amount: m.append("{}".format(self.Ceilometer_1_amount))
        else: m.append("")
        if not self.Ceilometer_1_height is None: m.append("{} (m)".format(self.Ceilometer_1_height))
        else: m.append("")

        if not self.Ceilometer_2_amount: m.append("{}".format(self.Ceilometer_2_amount))
        else: m.append("")
        if not self.Ceilometer_2_height is None: m.append("{} (m)".format(self.Ceilometer_2_height))
        else: m.append("")

        if not self.Ceilometer_3_amount: m.append("{}".format(self.Ceilometer_3_amount))
        else: m.append("")
        if not self.Ceilometer_3_height is None: m.append("{} (m)".format(self.Ceilometer_3_height))
        else: m.append("")

        if not self.Vis is None: m.append("{} (km)".format(self.Vis))
        else: m.append("")

        if not self.Pressure_MSLP is None: m.append("{} (hPa)".format(self.Pressure_MSLP))
        else: m.append("")
        if not self.Pressure_Station is None: m.append("{} (hPa)".format(self.Pressure_Station))
        else: m.append("")
        if not self.Pressure_QNH is None: m.append("{} (hPa)".format(self.Pressure_QNH))
        else: m.append("")

        return ",".join(m)

class BOM_AWS:
    '''
    Dictionary of AWS observations for a station
    '''
    def __init__(self, displayname: str, longitude: float, latitude: float):
        self.displayname = displayname
        self.longitude = longitude
        self.latitude = latitude
        self.filenameCSV:str = None
        self.observations: dict[datetime, AWS_Observation] = dict()
        self.xy = Point(longitude, latitude)
        self.label_latitude = latitude + 0.01
        self.label_longitude = longitude + 0.04

    def __str__(self):
        m = []
        m.append("BOM AWS: {} Lon:{} Lat: {} ".format(self.displayname, self.longitude, self.latitude))
        if not self.filenameCSV is None and self.filenameCSV: m.append("File: {}".format(self.filenameCSV))
        m.append("Observations count: {}".format(len(self.observations)))
        return "\n".join(m)
    
    def filter_nights_with_calm_winds(self, startdate: datetime, enddate: datetime, min_count) -> list[tuple[datetime, int, int]]:
        '''
        Filter nights with at least two observations of calm
        '''
        dt = startdate
        dt_last = dt
        count_all = 0
        count_calm = 0
        result: list[tuple[datetime, int, int]] = [] 
        while dt < enddate:
            if dt != dt_last and dt.hour == 0 and dt.minute == 0 and min_count <= count_calm :
                if not logging is None: logging.info("{}-{}-{}, {}, {}".format(dt_last.year, dt_last.month, dt_last.day, count_all, count_calm))
                result.append((datetime(dt_last.year, dt_last.month, dt_last.day, tzinfo=timezone.utc), count_all, count_calm))

            if dt.hour == 0 and dt.minute == 0:
                count_all = 0
                count_calm = 0
            
            if dt in self.observations.keys() and 9 <= dt.hour and dt.hour <= 20:
                obs = self.observations[dt]
                if not obs.Wind_Speed is None: 
                    count_all += 1
                    if obs.Wind_Speed == 0: 
                        count_calm += 1

            dt_last = dt
            dt += timedelta(minutes=30)
        
        if dt != dt_last and dt.hour == 0 and dt.minute == 0 and min_count <= count_calm :
            if not logging is None: logging.info(dt_last.year, dt_last.month, dt_last.day, count_all, count_calm)
            result.append((datetime(dt_last.year, dt_last.month, dt_last.day, tzinfo=timezone.utc), count_all, count_calm))

        return result

    def filter_relhum(self, startdate: datetime, enddate: datetime, rel_hum_min:int, min_count:int) -> list[tuple[datetime, int, int]]:
        '''
        Filter nights with at least two observations of calm
        '''
        dt = startdate
        dt_last = dt
        count_all = 0
        count_low_rh = 0
        result: list[tuple[datetime, int, int]] = [] 
        while dt < enddate:
            if dt != dt_last and dt.hour == 0 and dt.minute == 0 and min_count <= count_low_rh :
                if not logging is None: logging.info("{}-{}-{}, {}, {}".format(dt_last.year, dt_last.month, dt_last.day, count_all, count_low_rh))
                result.append((datetime(dt_last.year, dt_last.month, dt_last.day, tzinfo=timezone.utc), count_all, count_low_rh))

            if dt.hour == 0 and dt.minute == 0:
                count_all = 0
                count_low_rh = 0
            
            if dt in self.observations.keys() and 9 <= dt.hour and dt.hour <= 20:
                obs = self.observations[dt]
                if not obs.RelHum is None: 
                    count_all += 1
                    if obs.RelHum < rel_hum_min: 
                        count_low_rh += 1
                        if not logging is None: logging.info(obs)

            dt_last = dt
            dt += timedelta(minutes=30)
        
        if dt != dt_last and dt.hour == 0 and dt.minute == 0 and min_count <= count_low_rh :
            if not logging is None: logging.info("{}-{}-{}, {}, {}".format(dt_last.year, dt_last.month, dt_last.day, count_all, count_low_rh))
            result.append((datetime(dt_last.year, dt_last.month, dt_last.day, tzinfo=timezone.utc), count_all, count_low_rh))

        return result

    def get_obs(self, dt: datetime) -> AWS_Observation:
        '''
        Get observation

        Parameters
        ----------
        dt: datetime
            Date time of observation

        Returns
        -------
        obs: AWS_Observation or None
        '''
        if dt in self.observations.keys(): 
            return self.observations[dt]
        else: 
            return None

    def get_rose_bins(self):
        '''
        Determine the relative percentage of observation in each speed and direction bin
        '''
        total_count = 0
        # Initialized rose rose[bin, array[direction]]
        rose = []
        for t in spd_bins:
            rose.append(zeros(24))

        for o in self.observations.values():
            dir = o.Wind_Direction_Deg 
            if dir is None: dir = 0
            if dir == 360: dir = 0
            speed = o.Wind_Speed
            if speed is None: continue

            total_count += 1
            if speed == 0:
                for t in range(24): rose[0][t] += 1
            else:
                # Determine direction index
                d = 0
                if dir < 352.5:
                    for t in range (24):
                        if dir < 7.5 + t * 15: 
                            d = t
                            break

                # Determine speed bin index
                s = 0
                for t in range(1,len(spd_bins)):
                    if speed < spd_bins[t]:
                        s = t
                        break                    
                
                rose[s][d] += 1

        for s in range(len(spd_bins)):
            for d in range(24):
                 rose[s][d] = rose[s][d] / total_count 

        return rose

    def list_BOM_Moranbabah(self, date_start: datetime, date_end: datetime) -> list[str]:
        '''
        Return a CSV list of data for specified period
        Date, Air T C(deg), Dew Point C(deg), Wind Speed (km/h), Wind direction in degrees true
        '''
        arr = ["Date, Air T C(deg), Dew Point C(deg), Wind Speed (km/h), Wind direction in degrees true"]
        for dt in self.observations.keys():
            if dt < date_start: continue
            if date_end < dt: continue
            obs:AWS_Observation = self.observations[dt]
            arr.append("{}/{}/{} {}:{}, {}, {}, {}, {}".format(
                str(dt.day).zfill(2), str(dt.month).zfill(2), str(dt.year).zfill(4), str(dt.hour).zfill(2), str(dt.minute).zfill(2), 
                obs.Temperature_Air, obs.Temperature_DewPoint, obs.Wind_Speed, obs.Wind_Direction_Deg))
        return arr

    def plot_aws(self, ax: mpl.axes.Axes):
        '''
        Default way to plot a aws
        '''
        ax.plot(self.longitude, self.latitude, 'go', markersize=7)
        ax.text(self.label_longitude, self.label_latitude, self.displayname, color="black")
    
    def plot_wind_rose(self, filename: str):
        '''
        Define our wind rose function
        Step 8 Courtesy of https://gist.github.com/phobson/41b41bdd157a2bcf6e14
        '''
        rosedata = self.get_rose_bins()
        bar_width = 15 * pi / 180
        bar_dir = []
        for i in range(24): bar_dir.append(i * 15 * pi / 180)

        cmap = mpl.colormaps["tab10"]
        N = len(spd_bins)
        colors = cmap(linspace(0,1, N))
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
        fig.suptitle("{} wind rose 2018-01-01 and 2021-07-29".format(self.displayname))
        ax.set_theta_direction('clockwise')
        ax.set_theta_zero_location('N')
        ax.set_xticks([0,45 * pi / 180,90 * pi / 180,135 * pi / 180,180 * pi / 180,225 * pi / 180,270 * pi / 180,315 * pi / 180], ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])

        bottoms = []
        for n in range(N):
            if n == 0:
                # first column only
                bottoms = rosedata[n]
                ax.bar(bar_dir, bottoms, width=bar_width, color=colors[n], edgecolor='none', label="calm", linewidth=0)

            # all other columns
            else:
                l = "{} to {} (km/h)".format(spd_bins[n-1], spd_bins[n]) if n < N -1 else " above {} (km/h)".format(spd_bins[n-1])
                ax.bar(bar_dir, rosedata[n], width=bar_width, bottom=bottoms, color=colors[n], edgecolor='none', label=l, linewidth=0)
                for i in range(24): bottoms[i] += rosedata[n][i]

        ax.legend(loc=(0.8, 0.95), ncol=1)
       
        if not filename is None:
            plt.savefig(filename)
            m = "Chart generated: {}".format(filename)
            print(m)
            if not logger is None: logger.info(m)
        else:
            plt.show()

        plt.close("all")
        
class BOM_AWS_Moranbah(BOM_AWS):
    def __init__(self):
        super().__init__("Moranbah AWS", 148.075251, -22.062148)

def BOM_AWS_from_csv(displayname: str, filenameCSV: str) -> BOM_AWS:
    '''
    Factory method returning populated BOM object or None
    '''
    if not path.exists(filenameCSV):
        print("Error: File {} does not exist".format(filenameCSV))
        return None
    
    if displayname == "Moranbah AWS":
        bom = BOM_AWS_Moranbah()
        bom.filenameCSV = filenameCSV
        with open(filenameCSV, "r", newline="") as csvfile:
            fieldnames = [
                    "UTCDate",
                    "UTCTime",
                    "Precip10minutesmm",
                    "Precipsince9ammm",
                    "AirC",
                    "WetBulbC",
                    "DewPointC",
                    "RelHum%",
                    "WindSpeedkmh",
                    "WindDirectionDeg",
                    "WindGust",
                    "Ceilometer cloud amount (of first group)",
                    "Ceilometer cloud height (of first group) in meter",
                    "Ceilometer cloud amount (of second group)",
                    "Ceilometer cloud height (of second group) in meter",
                    "Ceilometer cloud amount (of third group)",
                    "Ceilometer cloud height (of third group) in meter",
                    "Viskm",
                    "MSLPhPa",
                    "StationPressurehPa",
                    "QNHhPa"
                ]
            reader = csv.DictReader(csvfile, fieldnames=fieldnames)
            i = 0
            for row in reader:
                if i > 0:
                    
                    utc_date = row["UTCDate"].strip()
                    if utc_date is None or not utc_date: continue
                    dateparts = utc_date.split("/")

                    utc_time = row["UTCTime"].strip()
                    if utc_time is None or not utc_time: continue
                    timeparts = utc_time.split(":")

                    date = datetime(int(dateparts[2]),int(dateparts[1]), int(dateparts[0]), int(timeparts[0]), int(timeparts[1]), 0, tzinfo=timezone.utc)
                    obs = AWS_Observation(displayname, date, bom.longitude, bom.latitude)

                    w_d = row["WindDirectionDeg"].strip()
                    if not w_d is None and w_d: obs.Wind_Direction_Deg = int(w_d)

                    w_g = row["WindGust"].strip()
                    if not w_g is None and w_g: obs.Wind_Gust = float(w_g)

                    w_s = row["WindSpeedkmh"].strip()
                    if not w_s is None and w_s: obs.Wind_Speed = float(w_s)

                    p10 = row["Precip10minutesmm"].strip()
                    if not p10 is None and p10: obs.Precip_10minutes = float(p10)

                    p9am = row["Precipsince9ammm"].strip()
                    if not p9am is None and p9am: obs.Precip_since9am = float(p9am)

                    t_a = row["AirC"].strip()
                    if not t_a is None and t_a:  obs.Temperature_Air = float(t_a)
                    
                    t_w = row["WetBulbC"].strip()
                    if not t_w is None and t_w: obs.Temperature_WetBulb = float(t_w)
                    
                    t_d = row["DewPointC"].strip()
                    if not t_d is None and t_d: obs.Temperature_DewPoint = float(t_d)

                    rh = row["RelHum%"].strip()
                    if not rh is None and rh: obs.RelHum = int(rh)
                    
                    obs.Ceilometer_1_amount = row["Ceilometer cloud amount (of first group)"].strip()

                    h1 = row["Ceilometer cloud height (of first group) in meter"].strip()
                    if not h1 is None and h1: obs.Ceilometer_1_height = int(h1)

                    obs.Ceilometer_2_amount = row["Ceilometer cloud amount (of second group)"].strip()

                    h2 = row["Ceilometer cloud height (of second group) in meter"].strip()
                    if not h2 is None and h2: obs.Ceilometer_2_height = int(h2)

                    obs.Ceilometer_3_amount = row["Ceilometer cloud amount (of third group)"].strip()

                    h3 = row["Ceilometer cloud height (of third group) in meter"].strip()
                    if not h3 is None and h3: obs.Ceilometer_3_height = int(h3)

                    v = row["Viskm"].strip()
                    if not v is None and v: obs.Vis = float(row["Viskm"].strip())

                    p_m = row["MSLPhPa"].strip()
                    if not p_m is None and p_m: obs.Pressure_MSLP = float(p_m)

                    p_s = row["StationPressurehPa"].strip()
                    if not p_s is None and p_s: obs.Pressure_Station = float(p_s)

                    p_q = row["QNHhPa"].strip()
                    if not p_q is None and p_q: obs.Pressure_QNH = float(p_q)
                    
                    bom.observations[date] = obs

                i += 1
        return bom
    else:
        print("Error: Unknown AWS {}".format(displayname))

def create_aws(name:str):
    if name == "Moranbah": return BOM_AWS_Moranbah()
    raise ValueError("Unknown name: {}".format(name))

def _handler_chart(args: Namespace):
    Config.create_log("BOM_AWS_Chart.log")
    bom = BOM_AWS_from_csv("Moranbah AWS", path_Moranbah_AWS_CSV)   
    bom.plot_wind_rose(None)

def _handler_list_observation(args: Namespace):
    '''
    List observation
    '''
    Config.create_log("BOM_AWS_List.log")
    
    dt_start = datetime(int(args.date_start[0:4]), int(args.date_start[4:6]), int(args.date_start[6:8]), int(args.date_start[8:10]), int(args.date_start[10:12]), tzinfo=timezone.utc)
    dt_end = datetime(int(args.date_end[0:4]), int(args.date_end[4:6]), int(args.date_end[6:8]), int(args.date_end[8:10]), int(args.date_end[10:12]), tzinfo=timezone.utc)
    bom = BOM_AWS_from_csv("Moranbah AWS", path_Moranbah_AWS_CSV )
    if dt_start != dt_end:
        for dt in bom.observations.keys():
            if dt < dt_start: continue
            if dt_end < dt: continue
            obs = bom.observations[dt]
            print(obs.to_csv())
            if not logging is None: logging.info(obs.to_csv)
    else:
        if dt_start in bom.observations.keys(): print(bom.observations[dt_start])

def _handler_list_period(args: Namespace):
    '''
    List observation
    '''
    Config.create_log("BOM_AWS_List.log")
    
    dt_end = datetime(int(args.date[0:4]), int(args.date[4:6]), int(args.date[6:8]), int(args.date[8:10]), tzinfo=timezone.utc)
    dt_start = dt_end - timedelta(hours=args.hours)
    bom = BOM_AWS_from_csv("Moranbah AWS", path_Moranbah_AWS_CSV )
    if dt_start != dt_end:
        for dt in bom.observations.keys():
            if dt < dt_start: continue
            if dt_end < dt: continue
            obs = bom.observations[dt]
            print(obs.to_csv())
            if not logging is None: logging.info(obs.to_csv)
    else:
        if dt_start in bom.observations.keys(): print(bom.observations[dt_start])

if __name__ == '__main__':
    parser = ArgumentParser(prog="BOM_AWS")
    subparsers = parser.add_subparsers(help="subcommand help", required=True)

    parser_chart = subparsers.add_parser("chart", help="Chart windrose for all observations in the file")
    parser_chart.add_argument("type", choices=["windrose"])
    parser_chart.set_defaults(func=_handler_chart)

    parser_list_observation = subparsers.add_parser("list", help="Print observation for AWS between start and end date")
    parser_list_observation.add_argument("aws", type=str, choices=["Moranbah"], help="AWS name")
    parser_list_observation.add_argument("date_start", type=str, help="First date inclusive in the format YYYYMMDDHHMM e.g. 201909150430")
    parser_list_observation.add_argument("date_end", type=str, help="End date inclusive in the format YYYYMMDDHHMM e.g. 201909150430")
    parser_list_observation.set_defaults(func=_handler_list_observation)

    parser_list_period = subparsers.add_parser("list_period", help="Print observation for AWS for specific number of hours before date")
    parser_list_period.add_argument("aws", type=str, choices=["Moranbah"], help="AWS name")
    parser_list_period.add_argument("date", type=str, help="End date year in the format YYYYMMDDHH e.g. 2019090204")
    parser_list_period.add_argument("hours", type=int, help="Count of previous hours to display. e.g. 24")
    parser_list_period.set_defaults(func=_handler_list_period)

    args = parser.parse_args()
    
    if "date" in args:
        if not 12 == len(args.date):
            print("date: {} must be in the format YYYYMMDDHHMM".format(args.date))
            exit()

    if "date_start" in args:
        if not 12 == len(args.date_start):
            print("date_start: {} must be in the format YYYYMMDDHHMM".format(args.date_start))
            exit()

    if "date_end" in args:
        if not 12 == len(args.date_end):
            print("date_end: {} must be in the format YYYYMMDDHHMM".format(args.date_end))
            exit()

    args.func(args)
