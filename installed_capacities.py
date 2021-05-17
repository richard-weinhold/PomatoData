
import os
import requests
import pandas as pd
from pathlib import Path

os.chdir(r'C:\Users\riw\Documents\repositories\pomato_data')
from auxiliary import get_countries_regions_ffe

def anymod_installed_capacities(anymod_result_path):

    raw = pd.read_csv(anymod_result_path)
    raw.loc[raw.technology.str.contains("offshore"), "group"] +=  "offshore" 
    raw.loc[raw.technology.str.contains("onshore"), "group"] += "onshore" 
    raw = raw.drop(["id", "Unnamed: 10"], axis=1)
    
    installed_capacity  = (raw[(raw.variable.isin(["capaConv", "capaStSize","capaStIn"])) & \
                              # (raw.country == "DE") & \
                              (raw.timestep_superordinate_dispatch.str.contains("2030"))
                              ].groupby(["country", "group", "variable"]).sum()
                           ).reset_index().pivot(index=["country", "group"], 
                                                 columns="variable", 
                                                 values="value")
    return installed_capacity
if __name__ == "__main__": 
    
    
    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    anymod_result_path = wdir.joinpath("data_in/results_summary_8days_grid_202105061657.csv")
    
    country_data, nuts_data = get_countries_regions_ffe()    

    wind_potentials = pd.read_csv(wdir.joinpath('data_out/wind_potentials.csv'), index_col=0).set_index("name_short")
    wind_potentials = wind_potentials.rename(columns = {"value": "pontential"})
    pv_potentials = pd.read_csv(wdir.joinpath('data_out/pv_potentials.csv'), index_col=0).set_index("name_short")
    pv_potentials = pv_potentials.rename(columns = {"value": "pontential"})

    wind_availabilites = pd.read_csv(wdir.joinpath('data_out/wind_availabilites.csv'), index_col=0)
    wind_availabilites = wind_availabilites.rename(columns = {"value": "availability"})

    pv_availabilites = pd.read_csv(wdir.joinpath('data_out/pv_availabilites.csv'), index_col=0)
    pv_availabilites = pv_availabilites.rename(columns = {"value": "availability"})

    tmp_pv = wind_availabilites[["nuts_id", "value"]].groupby("nuts_id").mean()
    tmp_nuts = nuts_data[nuts_data.name_short.isin(tmp_pv.index)]
    tmp_nuts["value"] = tmp_pv.loc[tmp_nuts.name_short].values
    tmp_nuts.plot(column="value", legend=True)

    installed_capacities = anymod_installed_capacities(anymod_result_path)
    
    # Capacity is distributed based on the potential in MWh
    
    countries = ["DE", "BE", "FR", "LU", "NL"]

    wind_capacities = wind_potentials.copy()
    
    