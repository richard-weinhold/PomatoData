

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
    
    # for col in header.columns:
    #     if any([c in col for c in zones.index]) and not any([e in col for e in exclude]) and "load" in col:
    #         read_cols.append(col)
    
    demand = pd.read_csv(wdir.joinpath("data_in/demand/time_series_60min_singleindex.csv"),
                              header=0, usecols=read_cols)
    demand = demand.rename(columns=col_rename)
    demand.utc_timestamp = pd.to_datetime(demand.utc_timestamp).astype('datetime64[ns]')
    return demand
    
if __name__ == "__main__": 
    
    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    demand = get_demand(wdir)
    demand.to_csv(wdir.joinpath('data_out/demand/demand.csv'))

