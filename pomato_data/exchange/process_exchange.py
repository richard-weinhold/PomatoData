

import logging
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from pomato_data.auxiliary import get_countries_regions_ffe


def process_physical_crossborder_flow_entso_e(wdir, year):
    country_data, nuts_data = get_countries_regions_ffe()    
    usecols = ["DateTime", "OutAreaTypeCode", "OutAreaName", "OutMapCode", 
               "InAreaTypeCode", "InAreaName", "InMapCode", "FlowValue"]
    file_dir = wdir.joinpath("data_in/exchange/physical_crossborder_flow")
    files = [file for file in file_dir.glob("*.csv") if str(year) in str(file)]

    # Load Files    
    physical_flow_df = pd.DataFrame()
    for file in files:
        ### Load Raw Data
        physical_flow_raw = pd.read_csv(file, header=0, encoding="UTF-16",
                               sep="\t", usecols=usecols)
        
        condition_1 = (physical_flow_raw.OutAreaTypeCode == "CTY") & \
                        (physical_flow_raw.InAreaTypeCode == "CTY")
                        
        condition_1_1 = (physical_flow_raw.OutAreaTypeCode == "CTA") & \
                        (physical_flow_raw.InAreaTypeCode == "CTA")
                        
        condition_2 = physical_flow_raw.OutMapCode.isin(country_data.index) & \
                       physical_flow_raw.InMapCode.isin(country_data.index)
                    
        physical_flow_raw_cty = physical_flow_raw.loc[condition_1 & condition_2].copy()
        physical_flow_raw_cta = physical_flow_raw.loc[condition_1_1 & condition_2].copy()
        
        cols = ["DateTime", "OutAreaName", "OutMapCode", 
                "InAreaName", "InMapCode", "FlowValue"]
        physical_flow_raw = pd.merge(physical_flow_raw_cty[cols], physical_flow_raw_cta[cols], on=cols[:-1], how="outer", suffixes=("", "_cta"))
        condition = (physical_flow_raw.FlowValue.isna())&(physical_flow_raw.FlowValue_cta.notna())

        physical_flow_raw.loc[condition, "FlowValue"] = physical_flow_raw.loc[condition, "FlowValue_cta"]     
        physical_flow_raw = physical_flow_raw.drop("FlowValue_cta", axis=1)
        physical_flow_df = pd.concat([physical_flow_df, physical_flow_raw])
    
    physical_flow_df.sort_values("DateTime", inplace=True)
    physical_flow_df = physical_flow_df.groupby(['DateTime', 'OutAreaName', 'OutMapCode', 'InAreaName', 'InMapCode']).sum().reset_index()
    physical_flow_df = physical_flow_df[~(physical_flow_df.OutMapCode == physical_flow_df.InMapCode)]
    physical_flow_df.columns = ["utc_timestep", "from_zone_name", "from_zone", "to_zone_name", "to_zone", "value"]
    
    return physical_flow_df

def process_commercial_exchange_entso_e(wdir, year):
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
        condition = (exchange_raw.Capacity.isna())&(exchange_raw.Capacity_cta.notna())

        exchange_raw.loc[condition, "Capacity"] = exchange_raw.loc[condition, "Capacity_cta"]     
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
    import pomato_data

    wdir = Path(pomato_data.__path__[0]).parent 
    year = 2019
    physical_crossborder_flow = process_physical_crossborder_flow_entso_e(wdir, year)
    commercial_exchange = process_commercial_exchange_entso_e(wdir, year)
    
    physical_crossborder_flow.loc[(physical_crossborder_flow.from_zone == "DE") &\
                              (physical_crossborder_flow.to_zone == "LU"), "value" ].max()
        
    commercial_exchange.loc[(commercial_exchange.from_zone == "LU") &\
                            (commercial_exchange.to_zone == "DE"), "value" ].max()
        
    physical_crossborder_flow.to_csv(wdir.joinpath(f"data_out/exchange/physical_crossborder_flow_{year}.csv"))
    commercial_exchange.to_csv(wdir.joinpath(f"data_out/exchange/commercial_exchange_{year}.csv"))
                                     
    # physical_crossborder_flow                         


