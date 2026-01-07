from argparse import ArgumentParser, Namespace
from datetime import datetime
from ERA5 import ERA5
from matplotlib.axes import Axes
from shapely import box, Point

class Source:
    Sources = ["BaralabaNorth", "Blackwater", "BlairAthol", "Bluff", "Broadmeadow", "Byerwen", 
               "Callide", "CamebyDowns", "CapcoalOpen", "CapcoalUnderground", "Carborough", "Caval", "Centurion", "Clermont", "Collinsville", "Commodore", "Cook", "Coppabella", "Curragh",
               "Daunia", "Dawson", "Drake",
               "EnshamOpen", "EnshamUnderground",
               "Foxleigh",
               "GoonyellaOpen", "GoonyellaUnderground", "Grosvenor",
               "HailCreek", 
               "IsaacPlains",
               "Jax", "Jeebropilly", "Jellinbah",
               "Kestrel", "Kogan",
               "LakeVermont",
               "Meandu", "Meteor", "Middlemount", "Millennium", "Minerva", "Moorvale", "MoranbahNorth",
               "NewAcland", "Newlands", 
               "OakyCreek",
               "PeakDowns", "Poitrel",
               "Rolleston", 
               "Saraji", "Sonoma", "SouthWalker",
               "Yarrabee"]
    '''
    Common properties of sources of emissions

    activity: dictionary
        Key - year of emissions. Value - Activity for this year (e.g. raw coal )
    case_name: str
        Name of the folder to keep data (e.g. HailCreek -> giving folder Data/HailCreek)
    display_name: str
        Name of source to be displayed on charts
    longitude: float
        Longitude of the source
    latitude: float
        Latitiude of the source
    '''
    def __init__(self, case_name: str, display_name: str, longitude: float, latitude: float):
        self.activity = {}
        self.case_name = case_name
        self.display_name = display_name
        self.xy = Point(longitude, latitude)
        self.gridpoint = ERA5.calculate_ERA5_NearestGridPoint(self.xy)
        self.latitude = latitude 
        self.longitude = longitude 
        self.label_latitude = latitude + 0.01
        self.label_longitude = longitude + 0.04
        self.open_year = None

    def __str__(self) -> str:
        a = ["Code name: {} Display name: {}".format(self.case_name, self.display_name)]
        a.append("Longitude: {}, Latitude: {}".format(self.xy.x, self.xy.y))
        a.append("ERA5 grid point: Longitude: {}, Latitude: {}".format(self.gridpoint.x, self.gridpoint.y))
        a.append("Activity data:")
        for (k,v) in self.activity.items():
            a.append("\t{}: {}".format(k,v))
        return "\n".join(a) 

    def plot_source(self, ax: Axes):
        '''
        Default way to plot a source
        '''
        ax.plot(self.longitude, self.latitude, 'ko', markersize=7)
        ax.text(self.label_longitude, self.label_latitude, self.display_name, color="black")

    @staticmethod
    def plot_source_in_extent(dt: datetime, minlon:float , minlat:float, maxlon:float, maxlat:float, ax: Axes):
        '''
        Chart all sources within extent which have data for the FY covering the date
        
        :param dt: Datetime of the image
        :type dt: datetime
        :param minlon: Minimum longitude
        :type minlon: float
        :param minlat: Minimum latitude
        :type minlat: float
        :param maxlon: Maximum longitude
        :type maxlon: float
        :param maxlat: Maximum latitude
        :type maxlat: float
        :param ax: Axes to draw onto
        :type ax: Axes
        '''
        b = box(minlon, minlat, maxlon, maxlat)
        # Convert year to financial year e.g. FY2019 = Jul 2018 - Jun 2019
        year = dt.year
        month = dt.month
        if 6 < month: year += 1
        for s in Source.Sources:
            sc = create_source(s)
            if not year in sc.activity: continue
            if not b.contains(sc.xy): continue
            sc.plot_source(ax)

class Source_Mine_BaralabaNorth(Source):    
    '''
    Baralaba North Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("BaralabaNorth", "Baralaba North", 149.7919444, -24.13027778)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 1.475749
        self.activity[2024] = 1.57599
        self.activity[2023] = 1.55253
        self.activity[2022] = 1.570677
        self.activity[2021] = 1.237306
        self.activity[2020] = 0.927759
        self.activity[2019] = 0.816068
        self.activity[2018] = 0.013104
        self.activity[2016] = 0.469932
        self.activity[2015] = 0.7306
        self.activity[2014] = 0.699961
        self.activity[2013] = 0.736142
        self.activity[2012] = 1.407801
        self.activity[2011] = 0.176856
        self.activity[2010] = 0.526575

class Source_Mine_Blackwater(Source):    
    '''
    Blackwater Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Blackwater", "Blackwater", 148.88, -23.57)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 14.179661
        self.activity[2024] = 12.569179
        self.activity[2023] = 12.510615
        self.activity[2022] = 14.288329
        self.activity[2021] = 15.216273
        self.activity[2020] = 12.893485
        self.activity[2019] = 15.261759
        self.activity[2018] = 13.678878
        self.activity[2017] = 16.581448
        self.activity[2016] = 16.524952
        self.activity[2015] = 14.843923
        self.activity[2014] = 14.842473
        self.activity[2013] = 11.864726
        self.activity[2012] = 10.138183
        self.activity[2011] = 10.419869
        self.activity[2010] = 12.898272

class Source_Mine_BlairAthol(Source):    
    '''
    BlairAthol Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("BlairAthol", "Blair Athol", 147.5266667, -22.67138889)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 2.046428
        self.activity[2024] = 1.995018
        self.activity[2023] = 2.473123
        self.activity[2022] = 2.781663
        self.activity[2021] = 2.554579
        self.activity[2020] = 3.051972
        self.activity[2019] = 2.829318
        self.activity[2018] = 1.613958
        self.activity[2013] = 1.51943
        self.activity[2012] = 2.916876
        self.activity[2011] = 3.601402
        self.activity[2010] = 9.697668

class Source_Mine_Bluff(Source):    
    '''
    Bluff Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Bluff", "Bluff", 149.0786443, -23.58974663)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2024] = 0.345596
        self.activity[2023] = 0.477284
        self.activity[2021] = 0.222674
        self.activity[2020] = 0.833995
        self.activity[2019] = 0.244355

class Source_Mine_Broadmeadow(Source):    
    '''
    Broadmeadow Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Broadmeadow", "Broadmeadow", 147.9819727, -21.8025693)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 0.072842
        self.activity[2024] = 1.681385
        self.activity[2023] = 1.016688

class Source_Mine_Byerwen(Source):    
    '''
    Byerwen Coal Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Byerwen", "Byerwen", 147.88576941524298, -21.370020381040423), 
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 5.01595
        self.activity[2024] = 6.976767
        self.activity[2023] = 6.252491
        self.activity[2022] = 6.197225
        self.activity[2021] = 7.464168
        self.activity[2020] = 5.380656
        self.activity[2019] = 4.982422
        self.activity[2018] = 2.21716

class Source_Mine_Callide(Source):    
    '''
    Callide Coal Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Callide", "Callide", 150.6173, -24.3268)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 10.228106
        self.activity[2024] = 8.740215
        self.activity[2023] = 7.971285
        self.activity[2022] = 6.904985
        self.activity[2021] = 8.203566
        self.activity[2020] = 11.174655
        self.activity[2019] = 10.808867
        self.activity[2018] = 10.156142
        self.activity[2017] = 7.035557
        self.activity[2016] = 8.068659
        self.activity[2015] = 7.644939
        self.activity[2014] = 6.615974
        self.activity[2013] = 7.210073
        self.activity[2012] = 7.957905
        self.activity[2011] = 7.739364
        self.activity[2010] = 8.521304

class Source_Mine_CamebyDowns(Source):    
    '''
    Cameby Downs Coal Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("CamebyDowns", "Cameby Downs", 150.31, -26.98)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 3.802548
        self.activity[2024] = 3.659978
        self.activity[2023] = 3.14566
        self.activity[2022] = 2.860475
        self.activity[2021] = 2.927747
        self.activity[2020] = 2.694264
        self.activity[2019] = 2.943899
        self.activity[2018] = 2.601234
        self.activity[2017] = 2.626007
        self.activity[2016] = 2.296133
        self.activity[2015] = 1.927518
        self.activity[2014] = 1.761535
        self.activity[2013] = 2.137602
        self.activity[2012] = 1.905356
        self.activity[2011] = 0.663339

class Source_Mine_CapcoalOpen(Source):    
    '''
    Capcoal Open Coal Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("CapcoalOpen", "Capcoal Open", 148.580506, -22.990997)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 4.821038
        self.activity[2024] = 3.993347
        self.activity[2023] = 3.9559
        self.activity[2022] = 3.255383
        self.activity[2021] = 4.746059
        self.activity[2020] = 5.036083
        self.activity[2019] = 4.686095
        self.activity[2018] = 4.349667
        self.activity[2017] = 4.16918
        self.activity[2016] = 5.840633
        self.activity[2015] = 7.159842
        self.activity[2014] = 7.307955
        self.activity[2013] = 7.362873
        self.activity[2012] = 8.175305
        self.activity[2011] = 4.968362
        self.activity[2010] = 6.464411

class Source_Mine_CapcoalUnderground(Source):    
    '''
    Capcoal Underground Coal Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("CapcoalUnderground", "Capcoal Underground", 148.580506, -22.990997)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 5.43821
        self.activity[2024] = 5.040798
        self.activity[2023] = 4.61602
        self.activity[2022] = 5.147648
        self.activity[2021] = 7.658537
        self.activity[2020] = 5.925393
        self.activity[2019] = 7.12786
        self.activity[2018] = 7.609721
        self.activity[2017] = 7.404229
        self.activity[2016] = 10.241899
        self.activity[2015] = 7.181958
        self.activity[2014] = 7.340032
        self.activity[2013] = 4.452776
        self.activity[2012] = 4.52792
        self.activity[2011] = 3.665887
        self.activity[2010] = 4.04465

class Source_Mine_Carborough(Source):    
    '''
    Carborough Downs Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Carborough", "Carborough Downs", 148.2094, -21.9502)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 1.498829
        self.activity[2024] = 2.05485
        self.activity[2023] = 3.155513
        self.activity[2022] = 2.823402
        self.activity[2021] = 4.561695
        self.activity[2020] = 2.804339
        self.activity[2019] = 2.799118
        self.activity[2018] = 2.382783
        self.activity[2017] = 2.961494
        self.activity[2016] = 2.910847
        self.activity[2015] = 4.108925
        self.activity[2014] = 3.122643
        self.activity[2013] = 3.115673
        self.activity[2012] = 2.137684
        self.activity[2011] = 1.962601
        self.activity[2010] = 1.934341

class Source_Mine_Caval(Source):    
    '''
    Caval Ridge Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Caval", "Caval Ridge", 148.11385, -22.17942)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 13.713365
        self.activity[2024] = 13.001105
        self.activity[2023] = 14.216291
        self.activity[2022] = 14.893889
        self.activity[2021] = 13.892186
        self.activity[2020] = 15.109871
        self.activity[2019] = 12.883808
        self.activity[2018] = 14.746458
        self.activity[2017] = 12.36578
        self.activity[2016] = 13.162818
        self.activity[2015] = 10.659832
        self.activity[2014] = 3.374347

class Source_Mine_Centurion(Source):    
    '''
    Centurion Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Centurion", "Centurion", 147.96, -21.66)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 0.498548
        self.activity[2024] = 0.007094
        self.activity[2019] = 0.246378
        self.activity[2018] = 3.578289
        self.activity[2017] = 2.460215
        self.activity[2016] = 2.725874
        self.activity[2015] = 3.555989
        self.activity[2014] = 1.523941
        self.activity[2013] = 3.197219
        self.activity[2012] = 1.681474
        self.activity[2011] = 2.1692
        self.activity[2010] = 3.220558

class Source_Mine_Clermont(Source):    
    '''
    Clermont Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Clermont", "Clermont", 147.6195648, -22.7201075)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 11.44449
        self.activity[2024] = 11.511608
        self.activity[2023] = 9.681281
        self.activity[2022] = 9.719213
        self.activity[2021] = 8.362532
        self.activity[2020] = 11.878171
        self.activity[2019] = 11.787169
        self.activity[2018] = 11.861485
        self.activity[2017] = 11.57199
        self.activity[2016] = 13.396764
        self.activity[2015] = 12.21301
        self.activity[2014] = 12.566191
        self.activity[2013] = 11.414911
        self.activity[2012] = 6.224372
        self.activity[2011] = 6.541501
        self.activity[2010] = 0.719392

class Source_Mine_Collinsville(Source):    
    '''
    Collinsville Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Collinsville", "Collinsville", 147.845, -20.5525)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 2.748576
        self.activity[2024] = 4.142294
        self.activity[2023] = 3.376096
        self.activity[2022] = 3.582914
        self.activity[2021] = 2.739176
        self.activity[2020] = 4.170709
        self.activity[2019] = 4.936278
        self.activity[2018] = 3.515123
        self.activity[2017] = 1.708952
        self.activity[2016] = 1.647226
        self.activity[2015] = 4.727554
        self.activity[2014] = 1.741867
        self.activity[2013] = 5.041889
        self.activity[2012] = 4.86694
        self.activity[2011] = 3.537041
        self.activity[2010] = 5.549436

class Source_Mine_Commodore(Source):    
    '''
    Commodore Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Commodore", "Commodore", 151.2554422, -27.9343501)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 3.014581
        self.activity[2024] = 3.278
        self.activity[2023] = 3.390753
        self.activity[2022] = 3.507034
        self.activity[2021] = 3.538101
        self.activity[2020] = 3.83749
        self.activity[2019] = 3.4174
        self.activity[2018] = 3.810086
        self.activity[2017] = 3.580869
        self.activity[2016] = 3.559217
        self.activity[2015] = 3.478022
        self.activity[2014] = 3.833833
        self.activity[2013] = 3.750844
        self.activity[2012] = 2.92696
        self.activity[2011] = 3.560053
        self.activity[2010] = 3.487858

class Source_Mine_Cook(Source):    
    '''
    Cook Colliery specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Cook", "Cook Colliery", 148.8, -23.6)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 0.287922
        self.activity[2024] = 0.099009
        self.activity[2023] = 0.172543
        self.activity[2022] = 0.021013
        self.activity[2020] = 0.282177
        self.activity[2019] = 0.538361
        self.activity[2018] = 0.129951
        self.activity[2017] = 0.524549
        self.activity[2016] = 0.79695
        self.activity[2015] = 1.544576
        self.activity[2014] = 0.476027
        self.activity[2013] = 0.427284
        self.activity[2012] = 0.598788
        self.activity[2011] = 0.695219
        self.activity[2010] = 0.626111

class Source_Mine_Coppabella(Source):    
    '''
    Coppabella Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Coppabella", "Coppabella", 148.446868, -21.851618)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 2.843341
        self.activity[2024] = 3.808373
        self.activity[2023] = 3.543527
        self.activity[2022] = 3.91761
        self.activity[2021] = 3.188777
        self.activity[2020] = 4.20506
        self.activity[2019] = 3.611035
        self.activity[2018] = 4.260718
        self.activity[2017] = 4.249084
        self.activity[2016] = 3.844321
        self.activity[2015] = 4.610616
        self.activity[2014] = 5.073218
        self.activity[2013] = 5.321862
        self.activity[2012] = 3.311492
        self.activity[2011] = 3.117576
        self.activity[2010] = 4.429857

class Source_Mine_Curragh(Source):    
    '''
    Curragh Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Curragh", "Curragh", 148.84, -23.49)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 11.672994
        self.activity[2024] = 12.851376
        self.activity[2023] = 12.730426
        self.activity[2022] = 12.161036
        self.activity[2021] = 16.374991
        self.activity[2020] = 14.489382
        self.activity[2019] = 15.908151
        self.activity[2018] = 14.841574
        self.activity[2017] = 14.919729
        self.activity[2016] = 12.835398
        self.activity[2015] = 15.157786
        self.activity[2014] = 14.525229
        self.activity[2013] = 12.913806
        self.activity[2012] = 11.957456
        self.activity[2011] = 9.849941
        self.activity[2010] = 11.29756

class Source_Mine_Daunia(Source):    
    '''
    Daunia Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Daunia", "Daunia", 148.287, -22.075)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 5.848
        self.activity[2024] = 4.746873
        self.activity[2023] = 5.069133
        self.activity[2022] = 3.650831
        self.activity[2021] = 4.814067
        self.activity[2020] = 5.513453
        self.activity[2019] = 5.451899
        self.activity[2018] = 6.09367
        self.activity[2017] = 5.990593
        self.activity[2016] = 6.060224
        self.activity[2015] = 5.670173
        self.activity[2014] = 5.354395
        self.activity[2013] = 1.251003

class Source_Mine_Dawson(Source):    
    '''
    Dawson North And Central Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Dawson", "Dawson North And Central", 149.77, -24.59)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 8.650744
        self.activity[2024] = 8.786474
        self.activity[2023] = 9.215373
        self.activity[2022] = 6.087247
        self.activity[2021] = 8.979575
        self.activity[2020] = 10.555571
        self.activity[2019] = 9.925466
        self.activity[2018] = 8.659475
        self.activity[2017] = 12.127655
        self.activity[2016] = 11.172833
        self.activity[2015] = 11.399332
        self.activity[2014] = 10.934242
        self.activity[2013] = 11.203474

class Source_Mine_Drake(Source):
    '''
    Northern Hub - Drake Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Drake", "Northern Hub - Drake", 147.812727, -20.705189)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 4.92531
        self.activity[2024] = 5.192481
        self.activity[2023] = 4.967919
        self.activity[2022] = 5.312472
        self.activity[2021] = 4.72648
        self.activity[2020] = 4.247085
        self.activity[2019] = 4.932996
        self.activity[2018] = 6.289679
        self.activity[2017] = 5.199889
        self.activity[2016] = 3.316756
        self.activity[2015] = 0.897326

class Source_Mine_EnshamOpen(Source):
    '''
    Ensham Open Cut Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("EnshamOpen", "Ensham Open Cut", 148.497455, -23.454464)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2024] = 0.756535
        self.activity[2023] = 0.05904
        self.activity[2022] = 0.721954
        self.activity[2021] = 0.691626
        self.activity[2020] = 0.874803
        self.activity[2019] = 0.650873
        self.activity[2018] = 0.945378
        self.activity[2017] = 0.958221
        self.activity[2016] = 1.679386
        self.activity[2015] = 2.56003
        self.activity[2014] = 3.512919
        self.activity[2013] = 4.078527
        self.activity[2012] = 4.399
        self.activity[2011] = 4.660259
        self.activity[2010] = 7.538704

class Source_Mine_EnshamUnderground(Source):
    '''
    Ensham Underground Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("EnshamUnderground", "Ensham Underground", 148.497455, -23.454464)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 3.943947
        self.activity[2024] = 3.438613
        self.activity[2023] = 2.714122
        self.activity[2022] = 1.44524
        self.activity[2021] = 2.601336
        self.activity[2020] = 4.104454
        self.activity[2019] = 4.527653
        self.activity[2018] = 4.459686
        self.activity[2017] = 3.962464
        self.activity[2016] = 3.31712
        self.activity[2015] = 2.263302
        self.activity[2014] = 1.633304
        self.activity[2013] = 0.955683

class Source_Mine_Foxleigh(Source):
    '''
    Foxleigh Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Foxleigh", "Foxleigh", 148.803889, -22.997778)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 3.408058
        self.activity[2024] = 3.376218
        self.activity[2023] = 3.470607
        self.activity[2022] = 3.419789
        self.activity[2021] = 3.040692
        self.activity[2020] = 3.376775
        self.activity[2019] = 4.390007
        self.activity[2018] = 3.74315
        self.activity[2017] = 3.986321
        self.activity[2016] = 4.03104
        self.activity[2015] = 4.29081
        self.activity[2014] = 3.494756
        self.activity[2013] = 3.733505
        self.activity[2012] = 3.411992
        self.activity[2011] = 2.361341
        self.activity[2010] = 2.656252

class Source_Mine_GoonyellaOpen(Source):
    '''
    Goonyella Riverside Open Cut Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("GoonyellaOpen", "Goonyella Riverside And Broadmeadow", 147.9596994501292, -21.80358433490622)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 11.740355
        self.activity[2024] = 12.38967
        self.activity[2023] = 15.146307
        self.activity[2022] = 13.757088
        self.activity[2021] = 18.656673
        self.activity[2020] = 16.642173
        self.activity[2019] = 17.467007
        self.activity[2018] = 16.224617
        self.activity[2017] = 14.974655
        self.activity[2016] = 17.75071
        self.activity[2015] = 15.34855
        self.activity[2014] = 17.511895
        self.activity[2013] = 15.2186
        self.activity[2012] = 12.68893
        self.activity[2011] = 15.518242
        self.activity[2010] = 16.138992

class Source_Mine_GoonyellaUnderground(Source):
    '''
    Goonyella Riverside Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("GoonyellaUnderground", "Goonyella Riverside And Broadmeadow", 147.9596994501292, -21.80358433490622)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 4.96002
        self.activity[2024] = 5.786198
        self.activity[2023] = 7.056095
        self.activity[2022] = 8.395057
        self.activity[2021] = 6.843411
        self.activity[2020] = 6.542735
        self.activity[2019] = 5.265511
        self.activity[2018] = 5.135977
        self.activity[2017] = 5.534594
        self.activity[2016] = 6.528215
        self.activity[2015] = 5.85932
        self.activity[2014] = 4.084879
        self.activity[2013] = 3.871383
        self.activity[2012] = 3.278733
        self.activity[2011] = 1.459222
        self.activity[2010] = 3.837427

class Source_Mine_Grosvenor(Source):
    '''
    Grosvenor Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Grosvenor", "Grosvenor", 148.0106, -21.9744)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2024] = 6.941619
        self.activity[2023] = 5.376078
        self.activity[2022] = 1.992752
        self.activity[2021] = 0.005933
        self.activity[2020] = 5.096847
        self.activity[2019] = 6.199142
        self.activity[2018] = 6.017762
        self.activity[2017] = 3.819068
        self.activity[2016] = 1.183411
        self.activity[2015] = 0.224074

class Source_Mine_HailCreek(Source):    
    '''
    Hail Creek mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("HailCreek", "Hail Creek", 148.3631389, -21.490934)
    
        # "Hail Creek FactSheet 2024 - final.pdf" from Glencore website
        self.open_year = 2003
        
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 10.807512
        self.activity[2024] = 9.592586
        self.activity[2023] = 10.092999
        self.activity[2022] = 10.783931
        self.activity[2021] = 10.594483
        self.activity[2020] = 9.468696
        self.activity[2019] = 7.660496
        self.activity[2018] = 10.203613
        self.activity[2017] = 9.195752
        self.activity[2016] = 10.228814
        self.activity[2015] = 11.614548
        self.activity[2014] = 11.911568
        self.activity[2013] = 12.063330
        self.activity[2012] = 12.817071
        self.activity[2011] = 11.983327
        self.activity[2010] = 12.554146

class Source_Mine_IsaacPlains(Source):
    '''
    Isaac Plains Complex specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("IsaacPlains", "Isaac Plains Complex", 148.16, -22.09)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 3.498907
        self.activity[2024] = 3.562296
        self.activity[2023] = 4.205046
        self.activity[2022] = 3.382315
        self.activity[2021] = 2.600271
        self.activity[2020] = 3.020412
        self.activity[2019] = 2.928712
        self.activity[2018] = 1.677807
        self.activity[2017] = 1.73707
        self.activity[2016] = 0.344703
        self.activity[2015] = 1.67018
        self.activity[2014] = 2.629967
        self.activity[2013] = 2.717582
        self.activity[2012] = 3.345687
        self.activity[2011] = 2.212276
        self.activity[2010] = 2.849146

class Source_Mine_Jax(Source):
    '''
    Northern Hub - Jax Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Jax", "Northern Hub - Jax", 147.860166, -20.713736)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 0.703587
        self.activity[2024] = 0.2648
        self.activity[2023] = 0.180924
        self.activity[2025] = 2.166749
        self.activity[2024] = 2.587424
        self.activity[2023] = 2.234311
        self.activity[2022] = 2.150849
        self.activity[2021] = 1.778437
        self.activity[2020] = 1.881646
        self.activity[2019] = 1.162692
        self.activity[2014] = 0.110725
        self.activity[2013] = 0.49698

class Source_Mine_Jeebropilly(Source):
    '''
    Jeebropilly Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Jeebropilly", "Jeebropilly", 152.6592, -27.6526)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2020] = 0.487693
        self.activity[2019] = 1.262099
        self.activity[2018] = 1.211851
        self.activity[2017] = 1.226294
        self.activity[2016] = 1.226744
        self.activity[2015] = 1.163656
        self.activity[2014] = 1.50057
        self.activity[2013] = 1.565582
        self.activity[2012] = 1.475845
        self.activity[2011] = 1.407669
        self.activity[2010] = 1.515545

class Source_Mine_Jellinbah(Source):
    '''
    Jellinbah Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Jellinbah", "Jellinbah", 148.9366, -23.394)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 6.672883
        self.activity[2024] = 6.125257
        self.activity[2023] = 5.249274
        self.activity[2022] = 2.778166
        self.activity[2021] = 5.586954
        self.activity[2020] = 4.709559
        self.activity[2019] = 5.285186
        self.activity[2018] = 4.495773
        self.activity[2017] = 5.610171
        self.activity[2016] = 5.153641
        self.activity[2015] = 5.147167
        self.activity[2014] = 5.004716
        self.activity[2013] = 4.946308
        self.activity[2012] = 4.803075
        self.activity[2011] = 4.092197
        self.activity[2010] = 4.508299

class Source_Mine_Kestrel(Source):
    '''
    Kestrel Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Kestrel", "Kestrel", 148.36, -23.23)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 6.025365
        self.activity[2024] = 7.739428
        self.activity[2023] = 7.463895
        self.activity[2022] = 7.391618
        self.activity[2021] = 7.475643
        self.activity[2020] = 8.174914
        self.activity[2019] = 7.757922
        self.activity[2018] = 6.344984
        self.activity[2017] = 6.172288
        self.activity[2016] = 4.991471
        self.activity[2015] = 3.752375
        self.activity[2014] = 2.567672
        self.activity[2013] = 3.293671
        self.activity[2012] = 5.474899
        self.activity[2011] = 4.460103
        self.activity[2010] = 4.936431

class Source_Mine_Kogan(Source):
    '''
    Kogan Creek Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Kogan", "Kogan Creek", 150.78, -26.93)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 2.215656
        self.activity[2024] = 2.538874
        self.activity[2023] = 2.086943
        self.activity[2022] = 2.588952
        self.activity[2021] = 1.759296
        self.activity[2018] = 1.783758
        self.activity[2017] = 2.708293
        self.activity[2016] = 2.15121
        self.activity[2015] = 2.660646
        self.activity[2014] = 2.461869
        self.activity[2013] = 2.558923
        self.activity[2012] = 1.946665
        self.activity[2011] = 2.525512
        self.activity[2010] = 2.15036

class Source_Mine_LakeVermont(Source):
    '''
    Lake Vermont Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("LakeVermont", "Lake Vermont", 148.4081, -22.4723)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 10.593207
        self.activity[2024] = 11.236347
        self.activity[2023] = 10.603634
        self.activity[2022] = 11.156475
        self.activity[2021] = 9.498174
        self.activity[2020] = 8.382135
        self.activity[2019] = 10.720365
        self.activity[2018] = 11.567648
        self.activity[2017] = 10.517068
        self.activity[2016] = 11.274715
        self.activity[2015] = 10.391608
        self.activity[2014] = 9.536098
        self.activity[2013] = 5.787559
        self.activity[2012] = 4.41963
        self.activity[2011] = 4.00267
        self.activity[2010] = 4.837828

class Source_Mine_Meandu(Source):
    '''
    Meandu Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Meandu", "Meandu", 151.9, -26.81666667)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 7.092086
        self.activity[2024] = 7.243326
        self.activity[2023] = 7.778163
        self.activity[2022] = 6.840659
        self.activity[2021] = 6.232684
        self.activity[2020] = 7.151744
        self.activity[2019] = 7.961919
        self.activity[2018] = 8.132043
        self.activity[2017] = 7.814977
        self.activity[2016] = 5.874786
        self.activity[2015] = 4.800283
        self.activity[2014] = 4.948059
        self.activity[2013] = 5.283161
        self.activity[2012] = 5.947657
        self.activity[2011] = 6.557772
        self.activity[2010] = 6.659619

class Source_Mine_Meteor(Source):
    '''
    Meteor Downs South Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Meteor", "Meteor Downs South", 148.385912, -24.433867)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 1.409742
        self.activity[2024] = 1.240391
        self.activity[2023] = 1.506637
        self.activity[2022] = 1.313371
        self.activity[2021] = 1.429798
        self.activity[2020] = 0.592005
        self.activity[2019] = 0.454229
        self.activity[2018] = 0.07853

class Source_Mine_Middlemount(Source):
    '''
    Middlemount Coal Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Middlemount", "Middlemount", 148.6386, -22.8814)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 3.957985
        self.activity[2024] = 3.755354
        self.activity[2023] = 3.070446
        self.activity[2022] = 4.270343
        self.activity[2021] = 4.750764
        self.activity[2020] = 2.92156
        self.activity[2019] = 4.460388
        self.activity[2018] = 5.305452
        self.activity[2017] = 5.121606
        self.activity[2016] = 5.471761
        self.activity[2015] = 5.303681
        self.activity[2014] = 4.913938
        self.activity[2013] = 2.327261
        self.activity[2012] = 1.752339
        self.activity[2011] = 0.089937

class Source_Mine_Millennium(Source):
    '''
    Millennium Coal Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Millennium", "Millennium", 148.252578, -22.015744)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2024] = 0.524882
        self.activity[2023] = 0.309609
        self.activity[2022] = 0.46543
        self.activity[2020] = 0.400612
        self.activity[2019] = 1.350518
        self.activity[2018] = 3.584283
        self.activity[2017] = 3.216802
        self.activity[2016] = 4.596807
        self.activity[2015] = 4.496724
        self.activity[2014] = 4.28402
        self.activity[2013] = 4.101826
        self.activity[2012] = 2.565564
        self.activity[2011] = 2.048483
        self.activity[2010] = 1.387512

class Source_Mine_Minerva(Source):
    '''
    Minerva Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Minerva", "Minerva", 148.0603, -23.9095)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2022] = 0.042047
        self.activity[2021] = 0.468049
        self.activity[2020] = 0.729765
        self.activity[2019] = 1.526704
        self.activity[2018] = 1.829852
        self.activity[2017] = 2.076445
        self.activity[2016] = 2.151935
        self.activity[2015] = 2.359272
        self.activity[2014] = 2.656694
        self.activity[2013] = 2.732424
        self.activity[2012] = 2.731628
        self.activity[2011] = 2.499378
        self.activity[2010] = 2.783916

class Source_Mine_Moorvale(Source):
    '''
    Moorvale Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Moorvale", "Moorvale", 148.3691, -21.9708)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 2.2343
        self.activity[2024] = 3.663958
        self.activity[2023] = 3.689346
        self.activity[2022] = 1.789768
        self.activity[2021] = 2.04969
        self.activity[2020] = 2.166217
        self.activity[2019] = 3.118094
        self.activity[2018] = 3.577809
        self.activity[2017] = 2.787666
        self.activity[2016] = 2.911265
        self.activity[2015] = 4.171469
        self.activity[2014] = 4.015055
        self.activity[2013] = 3.561877
        self.activity[2012] = 3.716619
        self.activity[2011] = 3.250238
        self.activity[2010] = 3.658319

class Source_Mine_MoranbahNorth(Source):    
    '''
    MoaranbahNorth mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("MoranbahNorth", "Moranbah North", 147.948876, -21.744594)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 2.269233
        self.activity[2024] = 3.554822
        self.activity[2023] = 5.880889
        self.activity[2022] = 4.124200
        self.activity[2021] = 5.849695
        self.activity[2020] = 7.325862
        self.activity[2019] = 7.744618
        self.activity[2018] = 9.429130
        self.activity[2017] = 7.785700
        self.activity[2016] = 5.938807
        self.activity[2015] = 7.062973
        self.activity[2014] = 6.618071
        self.activity[2013] = 6.169404
        self.activity[2012] = 3.434443
        self.activity[2011] = 5.346132
        self.activity[2010] = 4.199348

class Source_Mine_NewAcland(Source):
    '''
    New Acland Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("NewAcland", "New Acland", 151.6991, -27.27)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 6.094
        self.activity[2024] = 1.936
        self.activity[2022] = 0.977876
        self.activity[2021] = 3.926227
        self.activity[2020] = 5.779731
        self.activity[2019] = 8.984441
        self.activity[2018] = 9.288694
        self.activity[2017] = 9.287031
        self.activity[2016] = 9.245564
        self.activity[2015] = 10.144139
        self.activity[2014] = 7.484799
        self.activity[2013] = 8.886234
        self.activity[2012] = 10.322544
        self.activity[2011] = 9.155037
        self.activity[2010] = 9.494579

class Source_Mine_Newlands(Source):
    '''
    Newlands Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Newlands", "Newlands", 151.6991, -27.27)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2023] = 0.385172
        self.activity[2022] = 3.987619
        self.activity[2021] = 7.165654
        self.activity[2020] = 5.619942
        self.activity[2019] = 5.630744
        self.activity[2018] = 6.328799
        self.activity[2017] = 6.335703
        self.activity[2016] = 7.7534
        self.activity[2015] = 5.326775
        self.activity[2014] = 6.569743
        self.activity[2013] = 6.470112
        self.activity[2012] = 5.903613
        self.activity[2011] = 5.536529
        self.activity[2010] = 6.064754

class Source_Mine_OakyCreek(Source):
    '''
    Oaky Creek Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("OakyCreek", "Oaky Creek", 148.51295, -23.06342)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 5.3459
        self.activity[2024] = 6.087553
        self.activity[2023] = 5.51071
        self.activity[2022] = 5.769407
        self.activity[2021] = 5.3489
        self.activity[2020] = 6.463815
        self.activity[2019] = 6.57151
        self.activity[2018] = 5.746766
        self.activity[2017] = 8.923348
        self.activity[2016] = 8.868075
        self.activity[2015] = 8.616467
        self.activity[2014] = 10.600978
        self.activity[2013] = 11.817425
        self.activity[2012] = 11.963776
        self.activity[2011] = 13.243876
        self.activity[2010] = 14.164956

class Source_Mine_PeakDowns(Source):
    '''
    Peak Downs Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("PeakDowns", "Peak Downs", 148.1797, -22.2549)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 17.77571
        self.activity[2024] = 15.755505
        self.activity[2023] = 20.183395
        self.activity[2022] = 18.543975
        self.activity[2021] = 22.430126
        self.activity[2020] = 21.207575
        self.activity[2019] = 20.016342
        self.activity[2018] = 21.753506
        self.activity[2017] = 21.516076
        self.activity[2016] = 17.942049
        self.activity[2015] = 18.72108
        self.activity[2014] = 16.020212
        self.activity[2013] = 16.260221
        self.activity[2012] = 11.142325
        self.activity[2011] = 13.18837
        self.activity[2010] = 14.78638

class Source_Mine_Poitrel(Source):
    '''
    Poitrel Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Poitrel", "Poitrel", 148.2551254, -22.06650421)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 7.498741
        self.activity[2024] = 6.652584
        self.activity[2023] = 6.470649
        self.activity[2022] = 5.421983
        self.activity[2021] = 5.270765
        self.activity[2020] = 5.772193
        self.activity[2019] = 5.444997
        self.activity[2018] = 5.22428
        self.activity[2017] = 4.361813
        self.activity[2016] = 4.224826
        self.activity[2015] = 4.665752
        self.activity[2014] = 4.354469
        self.activity[2013] = 4.177733
        self.activity[2012] = 3.800417
        self.activity[2011] = 3.680797
        self.activity[2010] = 4.515682

class Source_Mine_Rolleston(Source):
    '''
    Rolleston Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Rolleston", "Rolleston", 148.4323, -24.3949)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 12.045045
        self.activity[2024] = 12.385695
        self.activity[2023] = 12.931587
        self.activity[2022] = 13.01697
        self.activity[2021] = 12.018336
        self.activity[2020] = 15.040236
        self.activity[2019] = 15.077646
        self.activity[2018] = 15.361209
        self.activity[2017] = 13.420437
        self.activity[2016] = 12.86955
        self.activity[2015] = 10.711291
        self.activity[2014] = 10.711575
        self.activity[2013] = 9.825951
        self.activity[2012] = 9.128982
        self.activity[2011] = 6.229544
        self.activity[2010] = 7.089281

class Source_Mine_Saraji(Source):
    '''
    Saraji Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Saraji", "Saraji", 148.29, -22.374)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 13.104971
        self.activity[2024] = 11.243284
        self.activity[2023] = 13.140015
        self.activity[2022] = 13.995076
        self.activity[2021] = 13.15153
        self.activity[2020] = 14.606833
        self.activity[2019] = 14.051941
        self.activity[2018] = 15.232977
        self.activity[2017] = 12.987779
        self.activity[2016] = 13.671868
        self.activity[2015] = 13.656529
        self.activity[2014] = 14.115471
        self.activity[2013] = 11.417588
        self.activity[2012] = 9.15425
        self.activity[2011] = 8.638008
        self.activity[2010] = 10.512975

class Source_Mine_Sonoma(Source):
    '''
    Northern Hub - Sonoma Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Sonoma", "Northern Hub - Sonoma", 147.86, -20.62)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2021] = 0.832104
        self.activity[2020] = 0.830819
        self.activity[2019] = 1.448954
        self.activity[2018] = 0.838242
        self.activity[2017] = 1.699505
        self.activity[2016] = 3.81112
        self.activity[2015] = 4.814854
        self.activity[2014] = 6.189626
        self.activity[2013] = 4.074544
        self.activity[2012] = 3.875884
        self.activity[2011] = 5.126085
        self.activity[2010] = 5.998901

class Source_Mine_SouthWalker(Source):
    '''
    South Walker Creek Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("SouthWalker", "South Walker Creek", 148.466821, -21.788705)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 8.096194
        self.activity[2024] = 8.278823
        self.activity[2023] = 8.026313
        self.activity[2022] = 8.007194
        self.activity[2021] = 6.628546
        self.activity[2020] = 7.188059
        self.activity[2019] = 7.729459
        self.activity[2018] = 8.108267
        self.activity[2017] = 7.375485
        self.activity[2016] = 7.174153
        self.activity[2015] = 6.198254
        self.activity[2014] = 6.49738
        self.activity[2013] = 6.109798
        self.activity[2012] = 5.965406
        self.activity[2011] = 4.523884
        self.activity[2010] = 4.811648

class Source_Mine_Yarrabee(Source):
    '''
    Yarrabee Mine specifications
    '''
    def __init__(self):
        '''
        Define source data
        Longitude and Latitude taken from maps.google.com
        '''
        super().__init__("Yarrabee", "Yarrabee", 149.026371, -23.317678)
        # Queensland Government Open Data Portal, FY year
        # Source: https://www.data.qld.gov.au/dataset/coal-industry-review-statistical-tables/resource/bab54159-f38b-4e6f-8652-4b04bca29139
        self.activity[2025] = 3.498192
        self.activity[2024] = 2.257483
        self.activity[2023] = 2.548488
        self.activity[2022] = 2.95027
        self.activity[2021] = 3.023392
        self.activity[2020] = 3.708726
        self.activity[2019] = 3.241634
        self.activity[2018] = 3.315867
        self.activity[2017] = 3.641199
        self.activity[2016] = 3.347345
        self.activity[2015] = 3.578484
        self.activity[2014] = 3.941137
        self.activity[2013] = 3.461826
        self.activity[2012] = 2.928851
        self.activity[2011] = 2.37149
        self.activity[2010] = 4.190178

def create_source(name: str) -> Source:
    '''
    Factory method for creating a source

    Parameters:
    -----------
    casename: str
        case name for source

    Returns:
    --------
    source: Source
        Source or none
    '''
    if name == "BaralabaNorth": return Source_Mine_BaralabaNorth()
    if name == "Blackwater": return Source_Mine_Blackwater()
    if name == "BlairAthol": return Source_Mine_BlairAthol()
    if name == "Bluff": return Source_Mine_Bluff()
    if name == "Broadmeadow": return Source_Mine_Broadmeadow() 
    if name == "Byerwen": return Source_Mine_Byerwen()
    if name == "Callide": return Source_Mine_Callide()
    if name == "CamebyDowns": return Source_Mine_CamebyDowns()
    if name == "CapcoalOpen": return Source_Mine_CapcoalOpen()
    if name == "CapcoalUnderground": return Source_Mine_CapcoalUnderground()
    if name == "Carborough": return Source_Mine_Carborough()
    if name == "Caval": return Source_Mine_Caval()
    if name == "Centurion": return Source_Mine_Centurion()
    if name == "Clermont": return Source_Mine_Clermont()
    if name == "Collinsville": return Source_Mine_Collinsville()
    if name == "Commodore": return Source_Mine_Commodore()
    if name == "Cook": return Source_Mine_Cook()
    if name == "Coppabella": return Source_Mine_Coppabella()
    if name == "Curragh": return Source_Mine_Curragh()
    if name == "Daunia": return Source_Mine_Daunia()
    if name == "Dawson": return Source_Mine_Dawson()
    if name == "Drake": return Source_Mine_Drake()
    if name == "EnshamOpen": return Source_Mine_EnshamOpen()
    if name == "EnshamUnderground": return Source_Mine_EnshamUnderground()
    if name == "Foxleigh": return Source_Mine_Foxleigh()
    if name == "GoonyellaOpen": return Source_Mine_GoonyellaOpen()
    if name == "GoonyellaUnderground": return Source_Mine_GoonyellaUnderground()
    if name == "Grosvenor": return Source_Mine_Grosvenor()
    if name == "HailCreek": return Source_Mine_HailCreek()
    if name == "IsaacPlains": return Source_Mine_IsaacPlains()
    if name == "Jax": return Source_Mine_Jax()
    if name == "Jeebropilly": return Source_Mine_Jeebropilly()
    if name == "Jellinbah": return Source_Mine_Jellinbah()
    if name == "Kestrel": return Source_Mine_Kestrel()
    if name == "Kogan": return Source_Mine_Kogan()
    if name == "LakeVermont": return Source_Mine_LakeVermont()
    if name == "Meandu": return Source_Mine_Meandu()
    if name == "Meteor": return Source_Mine_Meteor()
    if name == "Middlemount": return Source_Mine_Middlemount()
    if name == "Millennium": return Source_Mine_Millennium()
    if name == "Minerva": return Source_Mine_Minerva()
    if name == "Moorvale": return Source_Mine_Moorvale()
    if name == "MoranbahNorth": return Source_Mine_MoranbahNorth()
    if name == "NewAcland": return Source_Mine_NewAcland()
    if name == "Newlands": return Source_Mine_Newlands()
    if name == "OakyCreek": return Source_Mine_OakyCreek()
    if name == "PeakDowns": return Source_Mine_PeakDowns()
    if name == "Poitrel": return Source_Mine_Poitrel()
    if name == "Rolleston": return Source_Mine_Rolleston()
    if name == "Saraji": return Source_Mine_Saraji()
    if name == "Sonoma": return Source_Mine_Sonoma()
    if name == "SouthWalker": return Source_Mine_SouthWalker()
    if name == "Yarrabee": return Source_Mine_Yarrabee()

    raise ValueError("Unknown name: {}".format(name))

def _handler_list(args: Namespace) -> None:
    source = create_source(args.source)
    print(source)

if __name__ == '__main__':
    '''
    Entry point into this module
    '''
    parser = ArgumentParser(prog="Source", description="Utility to list available sources")
    subparsers = parser.add_subparsers()

    parser_list = subparsers.add_parser("list", help="List source")
    parser_list.add_argument("source", type=str, choices=Source.Sources, help="Source of emissions")
    parser_list.set_defaults(func=_handler_list)

    args = parser.parse_args()
    args.func(args)  