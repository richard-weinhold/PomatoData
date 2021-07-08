
import geopandas as gpd
import pandas as pd
import numpy as np
import atlite
import xarray as xr
from pathlib import Path
import os
import shapely

import logging
logging.basicConfig(level=logging.INFO)

os.chdir(r'C:\Users\riw\Documents\repositories\pomato_data')
from pomato_data.auxiliary import get_countries_regions_ffe, get_eez_ffe
os.environ['NUMEXPR_MAX_THREADS'] = '16'

def read_in_opsd_data(filepath, tech='Onshore'):

    input_df = pd.read_csv(filepath,
                           usecols=['electrical_capacity','energy_source_level_2','technology','nuts_3_region','lon','lat',
                                    'federal_state','commissioning_date','decommissioning_date','voltage_level'],
                          dtype={'energy_source_level_2':'str','technology':'str','nuts_3_region':'str','federal_state':'str',
                                 'commissioning_date':'str','decommissioning_date':'str','voltage_level':'str'})
    
    
    tech_df = input_df.loc[input_df['technology']==tech].rename(columns={'electrical_capacity':'capacity'}).copy()
    #keep entries which are not decommissioned and which have an nuts_3_region entry
    # print(onshore_df['decommissioning_date'].unique())
    tech_filter_df = tech_df.loc[(tech_df['decommissioning_date'].isna())&(tech_df['nuts_3_region'].notnull())].copy()
    #MW -> kW
    tech_filter_df['capacity'] = tech_filter_df['capacity']*1000
    tech_filter_df['country'] = np.array([(list(nuts3)[0] + list(nuts3)[1]) for nuts3 in tech_filter_df['nuts_3_region'].values])
    tech_filter_df = tech_filter_df.rename(columns={'technology':'type','commissioning_date':'installation_date'})
    return tech_filter_df

def categorise_wind_turbines(wind_df):
    
    turbine_categories = [
        dict(name='Vestas_V25_200kW', up=400),
        dict(name='Vestas_V47_660kW', up=700),
        dict(name='Bonus_B1000_1000kW', up=1100),
        dict(name='Suzlon_S82_1.5_MW', up=1600),
        dict(name='Vestas_V66_1750kW', up=1900),
        dict(name='Vestas_V80_2MW_gridstreamer', up=2200),
        dict(name='Siemens_SWT_2300kW', up=2500),
        dict(name='Vestas_V90_3MW', up=5000000)
    ]
    wind_categories = pd.DataFrame(turbine_categories)
    wind_categories['capacity']= 0
    for c in range(len(turbine_categories)):
        if c==0:
            low = 0
        else:
            low = turbine_categories[c-1]['up']
        if c==(len(turbine_categories)-1):
            high = turbine_categories[c]['up']
        else:
            high = turbine_categories[c]['up']
        wind_categories.loc[c,'capacity'] = wind_df.loc[wind_df['capacity'].between(low,high+1,inclusive=False),'capacity'].sum()
    wind_categories['capacity_relative'] = wind_categories['capacity']/wind_categories['capacity'].sum()
    return wind_categories


def prepare_cutout(weather_year, countries, cache_file_path, cache_file_name):
    # weather_year = '2019'
    # cache_file_path = wdir.joinpath("data_temp")
    # cache_file_name = "core"
    
    #Geometry shape for ERA5 cutout
    country_data, nuts_data = get_countries_regions_ffe()    
    
    # Dimensions of cutout (only relevant if it does not exist yes)
    x1, y1, x2, y2 = shapely.ops.cascaded_union(country_data.loc[countries, "geometry"].values).bounds
    # Path to store cutout: cutout_stor_path = 'data_temp\\cache_file_name
    cutout_stor_path = cache_file_path.joinpath(cache_file_name + '-' + str(weather_year))
           
    # Define cutout
    cutout = atlite.Cutout(path=str(cutout_stor_path),
                           module='era5',
                           x=slice(x1-.2, x2+.2), y=slice(y1-.2, y2+.2),
                           chunks={'time':100},
                           time=weather_year)
    cutout.prepare()
    return cutout 

def get_availabilities_atlite(cutout, countries, opsd_filepath):

    #read in onshore power plant data from OPSD
    onshore_df = read_in_opsd_data(opsd_filepath)
    #categorise onshore power plants in turbine categories
    data_category = categorise_wind_turbines(onshore_df)
    country_data, nuts_data = get_countries_regions_ffe()    

    cols = ["utc_timestamp", "nuts_id", "value"]
    wind = pd.DataFrame(columns=cols)
    pv = pd.DataFrame(columns=cols)
    for country in countries:
        print("Creating availability TS for ", country)
        tmp_nuts = nuts_data.loc[nuts_data.country == country, ["name_short", "geometry"]].set_index("name_short")
        wind_xarray = xr.Dataset()
        #Iterate over turbine categories
        for ind, dc in data_category.iterrows():
            name = f"< {dc['up']} kW"
            # Get feed in time-series for turbine category from cutout
            wind_xarray[name] = cutout.wind(turbine=dc['name'], 
                                            shapes=tmp_nuts['geometry'], 
                                            per_unit=True)
            # Multiply availability with weight (relative capacity of turbine category compared to all turbine categories)
            wind_xarray[name] = wind_xarray[name]*dc['capacity_relative']
        #sum up timeseries of all turbine categories    
        wind_xarray['total'] = sum(wind_xarray[c] for c in wind_xarray)
        #convert to dataframe
        tmp_wind_df = wind_xarray['total'].to_pandas()
        tmp_wind_df.columns
        # if country == "DK":
        #     tmp_wind_df["DK032"] = tmp_wind_df["DK031"]
            
        tmp_wind_df = tmp_wind_df.stack().reset_index()
        tmp_wind_df.columns = cols
        wind = pd.concat([wind, tmp_wind_df], ignore_index=True)        
        
        #PV panels orientated towards north, south, west and east
        tmp_pv_df_0 = cutout.pv(panel="CSi", orientation={'slope': 30., 'azimuth': 0.},
                                shapes=tmp_nuts['geometry'], per_unit=True)
        tmp_pv_df_90 = cutout.pv(panel="CSi", orientation={'slope': 30., 'azimuth': 90.}, 
                                 shapes=tmp_nuts['geometry'], per_unit=True)
        tmp_pv_df_180 = cutout.pv(panel="CSi", orientation={'slope': 30., 'azimuth': 180.}, 
                                  shapes=tmp_nuts['geometry'], per_unit=True)
        tmp_pv_df_270 = cutout.pv(panel="CSi", orientation={'slope': 30., 'azimuth': 270.}, 
                                  shapes=tmp_nuts['geometry'], per_unit=True)
        tmp_pv_df = (1/4)*tmp_pv_df_0+(1/4)*tmp_pv_df_90+(1/4)*tmp_pv_df_180+(1/4)*tmp_pv_df_270
        tmp_pv_df = tmp_pv_df.to_pandas()
        tmp_pv_df = tmp_pv_df.stack().reset_index()
        tmp_pv_df.columns = cols
        pv = pd.concat([pv, tmp_pv_df], ignore_index=True)        
        
    return wind, pv

def offshore_eez_atlite(cutout):
    eez = get_eez_ffe()
    eez = eez.set_index("id_region")
    turbine_type = "Vestas_V164_7MW_offshore"
    wind_offshore_availability = cutout.wind(turbine_type, 
                                             shapes=eez['geometry'], 
                                             per_unit=True)
    wind_offshore_availability = wind_offshore_availability.to_pandas()
    wind_offshore_availability = wind_offshore_availability.stack().reset_index()
    wind_offshore_availability.columns = ["utc_timestamp", "id_region", "value"]   
    return wind_offshore_availability


# %%

if __name__ == "__main__":
    import pomato_data
    
    wdir = Path(pomato_data.__path__[0]).parent 
    opsd_filepath = Path(r"C:\Users\riw\Documents\repositories\pomato_2030\res_capacity\data\renewable_power_plants_DE.csv")
    countries = ["DE", "BE", "FR", "LU", "NL", "CH", "AT", "CZ", "DK", "PL", "SE", "ES", "PT", "UK", "NO", "IT"]
    # countries = ["NO"]
    
    # wind, pv = get_availabilities_atlite(str(2019), wdir.joinpath("data_temp"), "core", opsd_filepath, countries)
    # # Save Resulting Tables. 
    # wind.to_csv(wdir.joinpath('data_out/res_availability/wind_availability_2019.csv'))
    # pv.to_csv(wdir.joinpath('data_out/res_availability/pv_availability_2019.csv'))
    
    # offshore = offshore_eez_ffe(str(2019), wdir.joinpath("data_temp"), "core")
    # offshore.to_csv(wdir.joinpath('data_out/res_availability/offshore_availability_2019.csv'))

    # %%
    
    import matplotlib.pyplot as plt
    
    wind = pd.read_csv(wdir.joinpath('data_out/res_availability/wind_availability_2019.csv'))
    pv = pd.read_csv(wdir.joinpath('data_out/res_availability/pv_availability_2019.csv'))
    
    fig, ax = plt.subplots(1, 1)
    pv = pv[["nuts_id", "value"]].groupby("nuts_id").mean()
    country_data, nuts_data = get_countries_regions_ffe()    
    pv = pd.merge(pv, nuts_data, left_index=True, right_on="name_short")
    gpd.GeoDataFrame(pv, geometry="geometry").plot(column="value", legend=True, ax=ax)  
    ax.axes.get_xaxis().set_ticks([])    
    ax.axes.get_yaxis().set_ticks([])    
    fig.tight_layout()    
    
    fig, ax = plt.subplots(1, 1)
    wind = wind[["nuts_id", "value"]].groupby("nuts_id").mean()
    country_data, nuts_data = get_countries_regions_ffe()    
    wind = pd.merge(wind, nuts_data, left_index=True, right_on="name_short")
    gpd.GeoDataFrame(wind, geometry="geometry").plot(column="value", legend=True, ax=ax)  
    ax.axes.get_xaxis().set_ticks([])    
    ax.axes.get_yaxis().set_ticks([])    
    fig.tight_layout()    









