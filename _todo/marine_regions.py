# -*- coding: utf-8 -*-
"""
Created on Thu May 20 12:15:42 2021

@author: s792
"""

import geopandas as gpd
import pandas as pd

marine_regions = gpd.read_file('World_EEZ_v11_20191118/eez_v11.shp')
iso3166_alpha3 = ['ALB','BEL','BGR','BIH','CYP','DEU','DNK','ESP','EST','FIN',
                  'FRA','GBR','GRC','HRV','IRL','ISL','ITA','LTU','LVA','MLT',
                  'MNE','NLD','NOR','POL','PRT','ROU','SVN','SWE','TUR','UKR']
regions = marine_regions.loc[marine_regions['ISO_TER1'].isin(iso3166_alpha3)]
