
import os
import requests
import pandas as pd
from pathlib import Path

os.chdir(r'C:\Users\riw\Documents\repositories\pomato_data')
from pomato_data.auxiliary import get_countries_regions_ffe

def get_pontetials_ffe():
    country_data, nuts_data = get_countries_regions_ffe()    

    # Download Potentials 
    url = "http://opendata.ffe.de:3000/v_opendata?id_opendata=eq.50"
    result = requests.get(url)
    data = pd.DataFrame.from_dict(result.json()[0]["data"])
    
    data = data[['id_region', 'internal_id_1', 'value']]

    wind = pd.merge(nuts_data, data[data.internal_id_1 == 1],  on="id_region", how="left")
    wind = wind.drop(["id_region", "internal_id_1"], axis=1).fillna(.1) 
    # default has top be non-zero for the capacity regionalisation. 
    
    pv = pd.merge(nuts_data, data[data.internal_id_1.isin([75, 76])],  on="id_region", how="left")
    pv = pv.drop(["id_region", "internal_id_1"], axis=1).fillna(.1)
    pv = pv.groupby(["name", "name_short", "country"]).sum().reset_index()
    pv = pd.merge(pv, nuts_data[["name_short", "geometry"]], on="name_short")
    
    return wind, pv


if __name__ == "__main__": 
    
    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    wind_potentials, pv_potentials = get_pontetials_ffe()
    wind_potentials.to_csv(wdir.joinpath('data_out/res_potential/wind_potential.csv'))
    pv_potentials.to_csv(wdir.joinpath('data_out/res_potential/pv_potential.csv'))
