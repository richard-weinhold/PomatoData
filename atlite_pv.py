# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 13:21:43 2021

@author: s792
"""

#%% Packages
import geopandas as gpd
import pandas as pd
import atlite
import xarray as xr
import pgeocode
import scipy.sparse as sp
import numpy as np
from pathlib import Path
import zipfile
import os
from collections import OrderedDict
from matplotlib import pyplot as plt 
import logging
logging.basicConfig(level=logging.INFO)

os.chdir('C:\\Users\\s792\\Documents\\pomato_2030\\res_availability\\atlite')
from atlite_auxiliary import *

#%% Working bench
os.chdir('C:\\Users\\s792\\Documents\\pomato_2030\\res_availability')
#### Preparation ####
#Used weather year
weather_year = '2012'

# #read in onshore power plant data from OPSD
# opsd_df = read_in_opsd_data(tech='Photovoltaics')
# opsd_df = opsd_df.rename(columns={'technology':'type','commissioning_date':'installation_date'})

#Geoemtry shape for ERA5 cutout
iso_a2_country = ['AT','BE','FR','DE','LU','NL']
# iso_a2_country = ['AT']
# iso_a2_country = ['BE','FR','DE','LU','NL']

#Geometry shape for splitting time series
#NUTS3 regions for central western europe (CWE)
nuts3_cwe_bd = gpd.read_file("data_temp/nuts3_cwe_boundaries.geojson")
#NUTS0 boundaries
nuts0_cwe_bd = gpd.read_file("data_temp/nuts0_cwe_boundaries.geojson")

#### Loop ####
for cntr in iso_a2_country:

    # cntr_record = list(filter(lambda c: c.attributes['ISO_A2'] == cntr, shp.records()))[0]
    # cntr_gdf = gpd.GeoDataFrame({**cntr_record.attributes, 'geometry':cntr_record.geometry},index=[0])
    # x1, y1, x2, y2 = cntr_gdf.loc[0]['geometry'].bounds
    if cntr=='FR':
        #FR is MultiPolygon with oversea islands
        x1, y1, x2, y2 = nuts0_cwe_bd.loc[nuts0_cwe_bd['CNTR_CODE']==cntr]['geometry'].iloc[0][0].bounds
    else:
        x1, y1, x2, y2 = nuts0_cwe_bd.loc[nuts0_cwe_bd['CNTR_CODE']==cntr]['geometry'].bounds.values[0]
    #Path to store cutout
    cutout_stor_path = 'data_temp\\'+cntr+'-'+weather_year
    # Define cutout
    cutout = atlite.Cutout(cutout_stor_path,
                           module='era5',
                           x=slice(x1-.2,x2+.2), y=slice(y1-.2, y2+.2),
                           chunks={'time':100},
                           time=weather_year)
    # Prepare cutout
    cutout.prepare()
    
    #NUTS3 regions for country
    nuts3_cntr_bd = nuts3_cwe_bd.loc[nuts3_cwe_bd['CNTR_CODE']==cntr].set_index('NUTS_ID')
    if cntr=='FR':
        #drop oversea regions
        nuts3_cntr_bd.drop(nuts3_cntr_bd.loc[nuts3_cntr_bd.index.str.contains('FRY')].index, inplace=True)
    
    ts = cutout.pv(panel="CSi", orientation={'slope': 25., 'azimuth': 0.}, shapes= nuts3_cntr_bd['geometry'], per_unit=True)
    
    #convert to dataframe
    ts_df = ts.to_pandas()
    #Store time-series
    ts_df.to_csv('data_out/atlite_ts_pv_'+cntr+'_'+weather_year+'.csv')
