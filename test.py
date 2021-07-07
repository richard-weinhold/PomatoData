# -*- coding: utf-8 -*-
"""
Created on Tue May 25 15:21:26 2021

@author: s792
"""

import geopandas as gpd
import atlite
import shapely
from pomato import POMATO
import os
from pathlib import Path
# os.chdir(r"C:\Users\s792\Documents\pomato_data\pomato_data")
# import sys
# sys.path.append("C:\\Users\\s792\\Documents\\pomato_data\\pomato_data")
# from auxiliary import get_countries_regions_ffe

#%%
# os.chdir(r"C:\Users\s792\Documents\pomato_data")
# hybas = gpd.read_file("data_in/res/hybas_eu/hybas_lake_eu_lev12_v1c/hybas_lake_eu_lev12_v1c.shp")
# hybas['geometry'].plot()

# hybas_data = hybas.drop(columns=['geometry'])

from pathlib import Path
import pandas as pd
import numpy as np
import requests 
import json 

def get_countries_regions_ffe(force_recalc=False):
    # Download the region types
    filepath = Path(__file__).parent.parent
    if filepath.joinpath("data_out/zones/zones.csv").is_file() and not force_recalc:
        zones = pd.read_csv(filepath.joinpath("data_out/zones/zones.csv"), index_col=0)
        zones['geometry'] = zones['geometry'].apply(shapely.wkt.loads)
        zones = gpd.GeoDataFrame(zones, geometry="geometry")

    else:
        print("Downloading zone data from FFE ")
        url = "http://opendata.ffe.de:3000/rpc/map_region_type?idregiontype=35"
        zones = gpd.read_file(url)
        zones = zones[['name', 'name_short', 'area_m2', "geometry"]].set_index("name_short")
    
    if filepath.joinpath("data_out/zones/nuts_data.csv").is_file() and not force_recalc:
        nuts_data = pd.read_csv(filepath.joinpath("data_out/zones/nuts_data.csv"), index_col=0)
        nuts_data['geometry'] = nuts_data['geometry'].apply(shapely.wkt.loads)
        nuts_data = gpd.GeoDataFrame(nuts_data, geometry="geometry")
    else:
        print("Downloading nuts data from FFE ")
        url = "http://opendata.ffe.de:3000/rpc/map_region_type?idregiontype=38"
        nuts_data = gpd.read_file(url)
        nuts_data = nuts_data[['id_region', 'name', 'name_short', 'area_m2', "geometry"]]
        geometry = []
        for geom in nuts_data['geometry']:
            if geom.is_valid:
                geometry.append(geom)
            else:
                geometry.append(geom.buffer(0))
                
        nuts_data = gpd.GeoDataFrame(nuts_data, geometry=geometry)
        nuts_data["country"] = nuts_data.name_short.str[:2]
        
        coastal = pd.read_csv(filepath.joinpath("data_in/regions/nuts_3_coastal.csv"), index_col=0)
        nuts_data["coastal"] = False
        nuts_data.loc[nuts_data.name_short.isin(coastal.index), "coastal"] = True
    return zones, nuts_data

#%% Hydro atlite func

def get_hydro_atlite(mato, 
                     weather_year='2017',
                     cache_file_path="C:\\Users\\s792\\Documents\\pomato_2030\\res_availability\\data_temp",
                     cache_file_name="core",
                     countries=["DE"]):

    # weather_year = '2020'
    # wdir = os.Path("C:\\Users\\s792\\Documents\\pomato_2030\\res_availability")
    # countries = ["NO"]

    # cache_file_path = wdir.joinpath("data_temp")
    # cache_file_name = "core"

        
    #Geoemtry shape for ERA5 cutout
    country_data, nuts_data = get_countries_regions_ffe()    
    
    x1, y1, x2, y2 = shapely.ops.cascaded_union(country_data.loc[countries, "geometry"].values).bounds
    #Path to store cutout
    # cutout_stor_path = 'data_temp\\'+cntr+'-'+weather_year
    cutout_stor_path = cache_file_path.joinpath(cache_file_name + '-' + str(weather_year))
           
    # Define cutout
    cutout = atlite.Cutout(path=str(cutout_stor_path),
                           module='era5',
                           x=slice(x1-.2, x2+.2), y=slice(y1-.2, y2+.2),
                           chunks={'time':100},
                           time=weather_year)
    cutout.prepare()
    
    # cols = ["utc_timestamp", "nuts_id", "value"]
    # wind = pd.DataFrame(columns=cols)
    # pv = pd.DataFrame(columns=cols)
    for zone in countries:
        # zone = "BE"
        
        mato.data.plants["zone"] = mato.data.nodes.loc[mato.data.plants.node, "zone"].values
        tmp_plants = mato.data.plants[(mato.data.plants.technology == "ror")&(mato.data.plants.zone == zone)]
        
        # eps = 0.1
        # t = tmp_plants[(tmp_plants.lat < 40 )]
        
        # geometry = [shapely.geometry.Point(xy) for xy in zip(tmp_plants.lon, tmp_plants.lat)]
        # tmp_plants = gpd.GeoDataFrame(tmp_plants, crs="EPSG:4326", geometry=geometry)
        
        basins = gpd.read_file(r"C:\Users\s792\Documents\pomato_data\data_in\res\hybas_eu\hybas_eu_lev12_v1c.zip")
        
        tmp = cutout.hydro(tmp_plants, basins, smooth=True)
        hydro_timeseries = tmp.to_pandas()
        hydro_timeseries = hydro_timeseries.T.fillna(0)
        hydro_timeseries.plot()
        
        # 1 m^3 = 1000 kg * 10 m 9.81 m/s^2 = 1000*10*9.81 Ws = (1000*10*9.81)/(3600 * 1e9) GWh
        
        hydro_timeseries_wh = hydro_timeseries*1000*10*9.81 / (3600*1e6)
        hydro_timeseries_wh.plot()
        
        # tmp_plants.loc[:, ["lat", "lon"]] =  53.060159, 8.865241 
        
#%%
wdir = Path(r"C:\Users\s792\Documents\pomato\pomato_wdir")
mato = POMATO(wdir=wdir, options_file="profiles/de_uniform_pricing.json")
# mato.load_data('data_input/pomato_dataset_five_zones/')
mato.load_data('data_input\\dataset_de_2019')

ts = get_hydro_atlite(mato)
