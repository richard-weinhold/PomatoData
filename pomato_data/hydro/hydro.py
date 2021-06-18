
import geopandas as gpd
import pandas as pd
import numpy as np
import atlite
import xarray as xr
from pathlib import Path
import os
import shapely
import matplotlib.pyplot as plt

import logging
logging.basicConfig(level=logging.INFO)

os.chdir(r'C:\Users\riw\Documents\repositories\pomato_data')
from pomato_data.auxiliary import get_countries_regions_ffe, get_eez_ffe
os.environ['NUMEXPR_MAX_THREADS'] = '16'

def get_hydro_atlite(weather_year, cache_file_path, cache_file_name, zones):

    # weather_year = '2020'
    # wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    # cache_file_path = wdir.joinpath("data_temp")
    # cache_file_name = "core"
    # zones = ['LU', 'ES', 'SE', 'AT', 'BE', 'CZ', 'DK', 'FR', 'DE', 'IT', 'NL', 'NO', 'PL', 'CH', 'UK']

    
    #Path to store cutout
    cutout_stor_path = cache_file_path.joinpath(cache_file_name + '-' + str(weather_year))
           
    # Define cutout
    cutout = atlite.Cutout(path=str(cutout_stor_path),
                           module='era5',
                           chunks={'time':100},
                           time=weather_year)
    cutout.prepare()

    plants = pd.read_csv(wdir.joinpath("data_in/hydro/jrc-hydro-power-plant-database.csv"), index_col=0)
    plants = plants[(plants.installed_capacity_MW > 0)&(plants.country_code.isin(zones))]
    
    inflow_plants = plants[plants["type"].isin(["HDAM", "HPHS"])]
    geometry = [shapely.geometry.Point(xy) for xy in zip(inflow_plants.lon, inflow_plants.lat)]
    inflow_plants = gpd.GeoDataFrame(inflow_plants, crs="EPSG:4326", geometry=geometry)

    hydro_timeseries = pd.DataFrame()
    basins = []
    for level in ["04", "05", "06", "07"]:
        basin = gpd.read_file(rf"C:\Users\riw\Documents\repositories\pomato_data\_todo\hydro\hybas_eu_lev{level}_v1c.zip")        
        geometry = []
        for geom in basin['geometry']:
            if geom.is_valid:
                geometry.append(geom)
            else:
                geometry.append(geom.buffer(0))
        basins.append(gpd.GeoDataFrame(basin, geometry=geometry).set_crs("EPSG:4326"))
    
    for zone in zones:
        tmp_plants = inflow_plants[(plants.country_code == zone)]
        if not tmp_plants.empty:
            tmp_plants_to_basin = gpd.sjoin(tmp_plants, basins[0], how='left', op='within')
            tmp_plants = tmp_plants[tmp_plants_to_basin.HYBAS_ID.notna()]
            
            inflow = pd.DataFrame()
            for basin in basins:
                tmp = cutout.hydro(tmp_plants, basin, smooth=True)
                if inflow.empty:
                    inflow = tmp.to_pandas().fillna(0)
                else:
                    inflow = inflow + tmp.to_pandas().fillna(0)
            hydro_timeseries = pd.concat([hydro_timeseries, inflow.T.fillna(0)], axis=1)

    cols = ["storage_capacity_MWh", "country_code", "installed_capacity_MW", "type", "avg_annual_generation_GWh"]
    flh_calc = inflow_plants.loc[inflow_plants.avg_annual_generation_GWh.notna(), cols].copy()
    flh_calc["flh"] = (1000*flh_calc.avg_annual_generation_GWh)/(flh_calc.installed_capacity_MW)
    flh_calc = flh_calc[["type", "flh", "country_code"]].groupby(["type", "country_code"]).mean()
    plants["flh"] = 0 
        
    for p in inflow_plants.iterrows():
        if (p[1].type, p[1].country_code) in flh_calc.index:
            inflow_plants.loc[p[0], "flh"] = flh_calc.loc[(p[1].type, p[1].country_code), "flh"]
        else:
            inflow_plants.loc[p[0], "flh"] = flh_calc.loc[p[1].type, "flh"].mean()
    
    # flh_calc.plot.scatter(x="installed_capacity_MW", y="flh")
    # flh_calc.plot.scatter(x="avg_annual_generation_GWh", y="storage_capacity_MWh")
    cond = inflow_plants.avg_annual_generation_GWh.isna()
    inflow_plants.loc[cond, "avg_annual_generation_GWh"] = (inflow_plants.loc[cond, "flh"] * inflow_plants.loc[cond, "installed_capacity_MW"])/1000
    inflows = hydro_timeseries.copy()    
    for p in inflows.columns:
        inflows.loc[:, p] = inflows.loc[:, p] * inflow_plants.loc[p, "avg_annual_generation_GWh"]*1000/inflows.loc[:, p].sum()
    inflows.index = inflows.index.rename("utc_timestamp")
    
    col_dict = {"name": "name", "installed_capacity_MW": "g_max", "pumping_MW": "d_max", "type": "technology", 
                "country_code": "zone", "storage_capacity_MWh": "storage_capacity", "lat": "lat", 
                "lon": "lon"}
    
    plants = plants[col_dict.keys()].rename(columns=col_dict)
    cond = plants.technology.isin(["HDAM", "HPHS"]) & (plants.storage_capacity.notna())
    avg_gen_storage_size = (plants.loc[cond, "storage_capacity"] / plants.loc[cond, "g_max"]).mean()
    cond = plants.technology.isin(["HDAM", "HPHS"]) & (plants.storage_capacity.isna())
    plants.loc[cond, "storage_capacity"] = avg_gen_storage_size  * plants.loc[cond, "g_max"]

    
    plants.technology.unique()
    plants["fuel"] = "hydro"
    plants.loc[:, "technology"] = plants.technology.replace({"HDAM": "reservoir", "HPHS": "psp", "HROR": "ror"})
    plants.loc[:, "plant_type"] = plants.technology.replace({"reservoir": "hydro_res", "psp": "hydro_psp", "ror": "hydro_ror"})
    plants.loc[:, "eta"] =  plants.technology.replace({"reservoir": 1, "psp": 0.75, "ror": 1})
    plants.loc[:, "mc_el"] =  plants.technology.replace({"reservoir": 7, "psp": 10, "ror": 1})
    plants.loc[:, "mc_heat"] =  0
    plants.loc[:, "h_max"] =  0
    plants.loc[:, "status"] =  "online"
    plants.loc[:, "commissioned"] =  1900
    plants[["chp", "heatarea", "city", "postcode", "street", "company"]] = None
    plants["node"] = None
    return plants , inflows
        
if __name__ == "__main__": 
    weather_year = '2019'
    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    cache_file_path = wdir.joinpath("data_temp")
    cache_file_name = "core"
    zones = ['LU', 'ES', 'SE', 'AT', 'BE', 'CZ', 'DK', 'FR', 'DE', 'IT', 'NL', 'NO', 'PL', 'CH', 'UK']
    plants, inflows = get_hydro_atlite(weather_year, cache_file_path, cache_file_name, zones)
    
    # plants.zone.unique()
    plants.to_csv(wdir.joinpath("data_out/hydro/plants.csv"))
    inflows.to_csv(wdir.joinpath(f"data_out/hydro/inflows_{weather_year}.csv"))
    
    
    
    
    