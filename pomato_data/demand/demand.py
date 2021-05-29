

import os
import requests
import pandas as pd
from pathlib import Path

os.chdir(r'C:\Users\riw\Documents\repositories\pomato_data')
from pomato_data.auxiliary import get_countries_regions_ffe

def get_demand(wdir):
    zones, nuts_data = get_countries_regions_ffe()    
    header = pd.read_csv(wdir.joinpath("data_in/demand/time_series_60min_singleindex.csv"), header=0, nrows=1)
    read_cols = ["utc_timestamp"]
    col_rename = {}
    for c in zones.index:
        if c + "_load_actual_entsoe_transparency" in header.columns:
            read_cols.append(c + "_load_actual_entsoe_transparency")
            col_rename[c + "_load_actual_entsoe_transparency"] = c
    read_cols.append("GB_GBN_load_actual_entsoe_transparency")
    col_rename["GB_GBN_load_actual_entsoe_transparency"] = "UK"
    
    demand = pd.read_csv(wdir.joinpath("data_in/demand/time_series_60min_singleindex.csv"),
                              header=0, usecols=read_cols)
    demand = demand.rename(columns=col_rename)
    demand.utc_timestamp = pd.to_datetime(demand.utc_timestamp).astype('datetime64[ns]')
    return demand
    
def get_demand_entso_e(wdir):
    country_data, nuts_data = get_countries_regions_ffe()    
    usecols = ["DateTime", "AreaTypeCode", "MapCode", "TotalLoadValue"]
    
    file_dir = wdir.joinpath("data_in\demand\ensto-e")
    files = [file for file in file_dir.glob("*.csv")]

    # Load Files    
    demand = pd.DataFrame()
    for file in files:
        ### Load Raw Data
        demand_raw = pd.read_csv(file, header=0, encoding="UTF-16",
                                 sep="\t", usecols=usecols)
        demand_raw["MapCode"] = demand_raw["MapCode"].replace("GB", "UK")
        demand_raw.DateTime = pd.to_datetime(demand_raw.DateTime).astype('datetime64[ns]')
        demand_raw['DateTime'] = demand_raw['DateTime'].apply(lambda x:x.replace(minute=0))
        demand_raw = demand_raw.groupby(usecols[:-1]).mean().reset_index()
        
        condition_cty = (demand_raw.AreaTypeCode == "CTY") 
        condition_cta = (demand_raw.AreaTypeCode == "CTA")
        condition_2 = demand_raw.MapCode.isin(country_data.index)
                    
        demand_cty = demand_raw.loc[condition_cty & condition_2].copy()
        demand_cta = demand_raw.loc[condition_cta & condition_2].copy()
        
        cols = ["DateTime", "MapCode", "TotalLoadValue"]
        demand_raw = pd.merge(demand_cty[cols], demand_cta[cols], on=cols[:-1], how="outer", suffixes=("", "_cta"))
        cond = (demand_raw.TotalLoadValue.isna())&(demand_raw.TotalLoadValue_cta.notna())

        demand_raw.loc[cond, "TotalLoadValue"] = demand_raw.loc[cond, "TotalLoadValue_cta"]     
        demand_raw = demand_raw.drop("TotalLoadValue_cta", axis=1)
        demand = pd.concat([demand, demand_raw])
        
    demand.sort_values("DateTime", inplace=True)

    # Fix missing values with the one an hour before
    demand = demand.pivot(index="DateTime", columns = "MapCode", values="TotalLoadValue")
    for zone in demand.columns:
        for i in demand.loc[demand[zone].isna(), zone].index:
               demand.iloc[demand.index.get_loc(i)][zone] = demand.iloc[demand.index.get_loc(i) - 1][zone]
    demand = demand.reset_index().rename(columns={"DateTime": "utc_timestamp"})
    
    return demand 

if __name__ == "__main__": 
    
    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    demand = get_demand_entso_e(wdir)
    demand.to_csv(wdir.joinpath('data_out/demand/demand.csv'))
    
