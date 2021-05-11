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
import scipy.sparse as sp
import numpy as np
from pathlib import Path
import zipfile
import os
from collections import OrderedDict
import logging
logging.basicConfig(level=logging.INFO)

os.chdir('C:\\Users\\s792\\Documents\\pomato_2030\\res_availability\\atlite')
from atlite_auxiliary import *

#%% Working bench
os.chdir('C:\\Users\\s792\\Documents\\pomato_2030\\res_availability')
#### Preparation ####
#Used weather year
weather_year = '2012'

turbine_categories = [
        dict(name='Vestas_V25_200kW', up=400),
        dict(name='Vestas_V47_660kW', up=700),
        dict(name='Bonus_B1000_1000kW', up=1100),
        dict(name='Suzlon_S82_1.5_MW', up=1600),
        dict(name='Vestas_V66_1750kW', up=1900),
        dict(name='Vestas_V80_2MW_gridstreamer', up=2200),
        dict(name='Siemens_SWT_2300kW', up=2500),
        dict(name='Vestas_V90_3MW', up=5000000)
    ]

#read in onshore power plant data from OPSD
onshore_df = read_in_opsd_data()
onshore_df = onshore_df.rename(columns={'technology':'type','commissioning_date':'installation_date'})
#categorise onshore power plants in turbine categories
data_categ = categorise_wind_turbines(onshore_df, turbine_categories)
data_categ['capacity_relative'] = data_categ['capacity']/data_categ['capacity'].sum()

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
    
    wind = xr.Dataset()
    #Iterate over turbine categories
    for ind, dc in data_categ.iterrows():
        name = f"< {dc['up']} kW"
        # Get feed in time-series for turbine category from cutout
        wind[name] = cutout.wind(turbine=dc['name'], shapes= nuts3_cntr_bd['geometry'], per_unit=True)
        # Multiply availability with weight (relative capacity of turbine category compared to all turbine categories)
        wind[name] = wind[name]*dc['capacity_relative']
    #sum up timeseries of all turbine categories    
    wind['total'] = sum(wind[c] for c in wind)
    
    #convert to dataframe
    wind_df = wind['total'].to_pandas()
    #Store time-series
    wind_df.to_csv('data_out/atlite_ts_wind_'+cntr+'_'+weather_year+'.csv')
