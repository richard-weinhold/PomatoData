

import geopandas as gpd
import pandas as pd
import numpy as np
import atlite
import xarray as xr
from pathlib import Path
import os
from shapely.ops import cascaded_union
import matplotlib.pyplot as plt

import logging
logging.basicConfig(level=logging.INFO)

os.chdir(r'C:\Users\riw\Documents\repositories\pomato_data')
from pomato_data.auxiliary import get_countries_regions_ffe

def process_physical_crossborder_flow(wdir):
    country_data, nuts_data = get_countries_regions_ffe()    
    usecols = ["DateTime", "OutAreaTypeCode", "OutAreaName", "OutMapCode", 
               "InAreaTypeCode", "InAreaName", "InMapCode", "FlowValue"]
    file_dir = wdir.joinpath("data_in/exchange/physical_crossborder_flow")
    files = [file for file in file_dir.glob("*.csv")]

    # Load Files    
    pcbf_df = pd.DataFrame(columns=usecols)
    for file in files:
        ### Load Raw Data
        cbpf_raw = pd.read_csv(file, header=0, encoding="UTF-16",
                               sep="\t", usecols=usecols)
        
        condition_1 = (cbpf_raw.OutAreaTypeCode == "CTY") & \
                        (cbpf_raw.InAreaTypeCode == "CTY")
                    
        condition_2 = cbpf_raw.OutMapCode.isin(country_data.index) & \
                       cbpf_raw.InMapCode.isin(country_data.index)
                    
        cbpf_raw = cbpf_raw.loc[condition_1 & condition_2]
        pcbf_df = pd.concat([pcbf_df, cbpf_raw])
    
    pcbf_df.sort_values("DateTime", inplace=True)
    pcbf_df = pcbf_df.drop(["OutAreaTypeCode", "InAreaTypeCode"], axis=1)
    pcbf_df = pcbf_df.groupby(['DateTime', 'OutAreaName', 'OutMapCode', 'InAreaName', 'InMapCode']).sum().reset_index()
    pcbf_df = pcbf_df[~(pcbf_df.OutMapCode == pcbf_df.InMapCode)]
    pcbf_df.columns = ["utc_timestep", "from_zone_name", "from_zone", "to_zone_name", "to_zone", "value"]
    return pcbf_df

# %%
if __name__ == "__main__":
    
    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    physical_crossborder_flow = process_physical_crossborder_flow(wdir)
    
    physical_crossborder_flow.loc[(physical_crossborder_flow.from_zone == "SE") &\
                              (physical_crossborder_flow.to_zone == "DE"), "value" ].max()
        
    # physical_crossborder_flow.to_csv(wdir.joinpath("data_out/exchange/physical_crossborder_flow.csv"))
                                     
    # physical_crossborder_flow                         


