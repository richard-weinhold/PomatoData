
import os
from pathlib import Path

import geopandas as gpd
import pandas as pd
import shapely
from pomato_data.auxiliary import get_countries_regions_ffe

os.environ['NUMEXPR_MAX_THREADS'] = '16'

def process_hydro_plants_with_atlite_inflows(cutout, zones, hydrobasins_path):

    plants = pd.read_csv(wdir.joinpath("data_in/hydro/jrc-hydro-power-plant-database.csv"), index_col=0)
    plants = plants[(plants.installed_capacity_MW > 0)&(plants.country_code.isin(zones))]
    
    inflow_plants = plants[plants["type"].isin(["HDAM", "HPHS"])]
    geometry = [shapely.geometry.Point(xy) for xy in zip(inflow_plants.lon, inflow_plants.lat)]
    inflow_plants = gpd.GeoDataFrame(inflow_plants, crs="EPSG:4326", geometry=geometry)

    hydro_timeseries = pd.DataFrame()
    basins = []
    for level in ["04", "05", "06", "07"]:
        basin = gpd.read_file(hydrobasins_path.joinpath(f"hybas_eu_lev{level}_v1c.zip"))        
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
    condition = inflow_plants.avg_annual_generation_GWh.isna()
    inflow_plants.loc[condition, "avg_annual_generation_GWh"] = (inflow_plants.loc[condition, "flh"] * inflow_plants.loc[condition, "installed_capacity_MW"])/1000
    inflows = hydro_timeseries.copy()    
    for p in inflows.columns:
        inflows.loc[:, p] = inflows.loc[:, p] * inflow_plants.loc[p, "avg_annual_generation_GWh"]*1000/inflows.loc[:, p].sum()
    inflows.index = inflows.index.rename("utc_timestamp")
    
    col_dict = {"name": "name", "installed_capacity_MW": "g_max", "pumping_MW": "d_max", "type": "technology", 
                "country_code": "zone", "storage_capacity_MWh": "storage_capacity", "lat": "lat", 
                "lon": "lon"}
    
    plants = plants[col_dict.keys()].rename(columns=col_dict)
    condition = plants.technology.isin(["HDAM", "HPHS"]) & (plants.storage_capacity.notna())
    avg_gen_storage_size = (plants.loc[condition, "storage_capacity"] / plants.loc[condition, "g_max"]).mean()
    condition = plants.technology.isin(["HDAM", "HPHS"]) & (plants.storage_capacity.isna())
    plants.loc[condition, "storage_capacity"] = avg_gen_storage_size  * plants.loc[condition, "g_max"]

    
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
        

def process_storage_level_entso_e(wdir, year):
    
    country_data, nuts_data = get_countries_regions_ffe()    
    usecols = ["DateTime", "AreaTypeCode", "MapCode", "StoredEnergy"]
    file_dir = wdir.joinpath("data_in/hydro/storage_level")
    files = [file for file in file_dir.glob("*.csv") if str(year) in str(file)]

    # Load Files    
    storage_level = pd.DataFrame()
    for file in files:
        ### Load Raw Data
        
        storage_level_raw = pd.read_csv(file, header=0, sep="\t", usecols=usecols)
        
        storage_level_raw["MapCode"] = storage_level_raw["MapCode"].replace("GB", "UK")

        storage_level_raw.DateTime = pd.to_datetime(storage_level_raw.DateTime).astype('datetime64[ns]')


        condition_cty = (storage_level_raw.AreaTypeCode == "CTY") 
        condition_cta = (storage_level_raw.AreaTypeCode == "CTA")
        condition_2 = storage_level_raw.MapCode.isin(country_data.index)
                    
        storage_level_cty = storage_level_raw.loc[condition_cty & condition_2].copy()
        storage_level_cta = storage_level_raw.loc[condition_cta & condition_2].copy()
        
        cols = ["DateTime", "MapCode", "StoredEnergy"]
        storage_level_raw = pd.merge(storage_level_cty[cols], storage_level_cta[cols], on=cols[:-1], how="outer", suffixes=("", "_cta"))
        condition = (storage_level_raw.StoredEnergy.isna())&(storage_level_raw.StoredEnergy_cta.notna())
        storage_level_raw.loc[condition, "StoredEnergy"] = storage_level_raw.loc[condition, "StoredEnergy_cta"]     
        storage_level_raw = storage_level_raw.drop("StoredEnergy_cta", axis=1)
        storage_level = pd.concat([storage_level, storage_level_raw])

    
    storage_level.sort_values("DateTime", inplace=True)
    
    storage_level.columns = ["utc_timestamp", "zone", "storage_level"]
    
    return storage_level

# %%
if __name__ == "__main__":
    import pomato_data
     
    weather_year = '2019'
    wdir = Path(pomato_data.__path__[0]).parent 
    cache_file_path = wdir.joinpath("data_temp")
    cache_file_name = "core"
    zones = ['LU', 'ES', 'SE', 'AT', 'BE', 'CZ', 'DK', 'FR', 'DE', 'IT', 'NL', 'NO', 'PL', 'CH', 'UK']
    plants, inflows = process_hydro_plants_with_atlite_inflows(weather_year, cache_file_path, cache_file_name, zones)
    storage_level = process_storage_level_entso_e(wdir, weather_year)
    
    # plants.zone.unique()
    plants.to_csv(wdir.joinpath("data_out/hydro/plants.csv"))
    inflows.to_csv(wdir.joinpath(f"data_out/hydro/inflows_{weather_year}.csv"))
    storage_level.to_csv(wdir.joinpath(f"data_out/hydro/storage_level_{weather_year}.csv"))
    
    
    
    