

import geopandas as gpd
import pandas as pd
import numpy as np
import atlite
import xarray as xr
from pathlib import Path
import os
from shapely.ops import cascaded_union
import matplotlib.pyplot as plt
import sys

import logging
logging.basicConfig(level=logging.INFO)

# os.chdir(r'C:\Users\riw\Documents\repositories\pomato_data')
file_dir = os.path.dirname(os.path.abspath(__file__))
package_dir = os.path.dirname(os.path.dirname(file_dir))
sys.path.append(package_dir)
from pomato_data.auxiliary import get_countries_regions_ffe

def process_physical_crossborder_flow(wdir, year):
    country_data, nuts_data = get_countries_regions_ffe()    
    usecols = ["DateTime", "OutAreaTypeCode", "OutAreaName", "OutMapCode", 
               "InAreaTypeCode", "InAreaName", "InMapCode", "FlowValue"]
    file_dir = wdir.joinpath("data_in/exchange/physical_crossborder_flow")
    files = [file for file in file_dir.glob("*.csv") if str(year) in str(file)]

    # Load Files    
    pcbf_df = pd.DataFrame()
    for file in files:
        ### Load Raw Data
        cbpf_raw = pd.read_csv(file, header=0, encoding="UTF-16",
                               sep="\t", usecols=usecols)
        
        condition_1 = (cbpf_raw.OutAreaTypeCode == "CTY") & \
                        (cbpf_raw.InAreaTypeCode == "CTY")
                        
        condition_1_1 = (cbpf_raw.OutAreaTypeCode == "CTA") & \
                        (cbpf_raw.InAreaTypeCode == "CTA")
                        
        condition_2 = cbpf_raw.OutMapCode.isin(country_data.index) & \
                       cbpf_raw.InMapCode.isin(country_data.index)
                    
        cbpf_raw_cty = cbpf_raw.loc[condition_1 & condition_2].copy()
        cbpf_raw_cta = cbpf_raw.loc[condition_1_1 & condition_2].copy()
        
        cols = ["DateTime", "OutAreaName", "OutMapCode", 
                "InAreaName", "InMapCode", "FlowValue"]
        cbpf_raw = pd.merge(cbpf_raw_cty[cols], cbpf_raw_cta[cols], on=cols[:-1], how="outer", suffixes=("", "_cta"))
        cond = (cbpf_raw.FlowValue.isna())&(cbpf_raw.FlowValue_cta.notna())

        cbpf_raw.loc[cond, "FlowValue"] = cbpf_raw.loc[cond, "FlowValue_cta"]     
        cbpf_raw = cbpf_raw.drop("FlowValue_cta", axis=1)
        pcbf_df = pd.concat([pcbf_df, cbpf_raw])
    
    pcbf_df.sort_values("DateTime", inplace=True)
    pcbf_df = pcbf_df.groupby(['DateTime', 'OutAreaName', 'OutMapCode', 'InAreaName', 'InMapCode']).sum().reset_index()
    pcbf_df = pcbf_df[~(pcbf_df.OutMapCode == pcbf_df.InMapCode)]
    pcbf_df.columns = ["utc_timestep", "from_zone_name", "from_zone", "to_zone_name", "to_zone", "value"]
    
    return pcbf_df

def process_commercial_exchange(wdir, year):
    country_data, nuts_data = get_countries_regions_ffe()    
    usecols = ["DateTime", "OutAreaTypeCode", "OutAreaName", "OutMapCode", 
               "InAreaTypeCode", "InAreaName", "InMapCode", "Capacity"]
    
    file_dir = wdir.joinpath("data_in/exchange/commercial_exchange")
    files = [file for file in file_dir.glob("*.csv") if str(year) in str(file)]

    # Load Files    
    exchange = pd.DataFrame()
    
    for file in files:
        1
        ### Load Raw Data
        exchange_raw = pd.read_csv(file, header=0, encoding="UTF-16",
                                   sep="\t", usecols=usecols)
        
        condition_1 = (exchange_raw.OutAreaTypeCode == "CTY") & \
                        (exchange_raw.InAreaTypeCode == "CTY")
        condition_1_1 = (exchange_raw.OutAreaTypeCode == "CTA") & \
                        (exchange_raw.InAreaTypeCode == "CTA")
                    
        condition_2 = exchange_raw.OutMapCode.isin(country_data.index) & \
                       exchange_raw.InMapCode.isin(country_data.index)
                    
        exchange_cty = exchange_raw.loc[condition_1 & condition_2].copy()
        exchange_cta = exchange_raw.loc[condition_1_1 & condition_2].copy()
        
        cols = ["DateTime", "OutAreaName", "OutMapCode", 
                "InAreaName", "InMapCode", "Capacity"]
        exchange_raw = pd.merge(exchange_cty[cols], exchange_cta[cols], on=cols[:-1], how="outer", suffixes=("", "_cta"))
        cond = (exchange_raw.Capacity.isna())&(exchange_raw.Capacity_cta.notna())

        exchange_raw.loc[cond, "Capacity"] = exchange_raw.loc[cond, "Capacity_cta"]     
        exchange_raw = exchange_raw.drop("Capacity_cta", axis=1)
        exchange = pd.concat([exchange, exchange_raw])
      
    exchange.sort_values("DateTime", inplace=True)
    # exchange = exchange.drop(["OutAreaTypeCode", "InAreaTypeCode"], axis=1)
    exchange = exchange.groupby(['DateTime', 'OutAreaName', 'OutMapCode', 'InAreaName', 'InMapCode']).sum().reset_index()
    exchange = exchange[~(exchange.OutMapCode == exchange.InMapCode)]
    exchange.columns = ["utc_timestep", "from_zone_name", "from_zone", "to_zone_name", "to_zone", "value"]
    
    return exchange

# %%
if __name__ == "__main__":
    
    wdir = Path(package_dir)
    year = 2020
    physical_crossborder_flow = process_physical_crossborder_flow(wdir, year)
    commercial_exchange = process_commercial_exchange(wdir, year)
    
    physical_crossborder_flow.loc[(physical_crossborder_flow.from_zone == "DE") &\
                              (physical_crossborder_flow.to_zone == "LU"), "value" ].max()
        
    commercial_exchange.loc[(commercial_exchange.from_zone == "LU") &\
                            (commercial_exchange.to_zone == "DE"), "value" ].max()
        
    physical_crossborder_flow.to_csv(wdir.joinpath(f"data_out/exchange/physical_crossborder_flow_{year}.csv"))
    commercial_exchange.to_csv(wdir.joinpath(f"data_out/exchange/commercial_exchange_{year}.csv"))
                                     
    # physical_crossborder_flow                         



# %%
