# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:19:53 2021

@author: s792
"""

#%% Packages
import pandas as pd
import geopandas as gpd
import pgeocode
import numpy as np
from collections import OrderedDict

#%% Data read-in and filtering
def read_in_opsd_data(path='C:\\Users\\s792\\Documents\\pomato_2030\\res_availability\\data_in\\renewable_power_plants_DE.csv',tech='Onshore'):
    input_df = pd.read_csv(path,
                           usecols=['electrical_capacity','energy_source_level_2','technology','nuts_3_region','lon','lat',
                                    'federal_state','commissioning_date','decommissioning_date','voltage_level'],
                          dtype={'energy_source_level_2':'str','technology':'str','nuts_3_region':'str','federal_state':'str',
                                 'commissioning_date':'str','decommissioning_date':'str','voltage_level':'str'})
    
    
    tech_df = input_df.loc[input_df['technology']==tech].rename(columns={'electrical_capacity':'capacity'}).copy()
    #keep entries which are not decommissioned and which have an nuts_3_region entry
    # print(onshore_df['decommissioning_date'].unique())
    tech_fltr_df = tech_df.loc[(tech_df['decommissioning_date'].isna())&(tech_df['nuts_3_region'].notnull())].copy()
    #MW -> kW
    tech_fltr_df['capacity'] = tech_fltr_df['capacity']*1000
    tech_fltr_df['country'] = np.array([(list(nuts3)[0] + list(nuts3)[1]) for nuts3 in tech_fltr_df['nuts_3_region'].values])
    return tech_fltr_df

# onshore_df = read_in_opsd_data()
#%% Grouping
#Wind turbine categories
# turbine_categories = [
#         dict(name='Vestas_V25_200kW', up=400),
#         dict(name='Vestas_V47_660kW', up=700),
#         dict(name='Bonus_B1000_1000kW', up=1100),
#         dict(name='Suzlon_S82_1.5_MW', up=1600),
#         dict(name='Vestas_V66_1750kW', up=1900),
#         dict(name='Vestas_V80_2MW_gridstreamer', up=2200),
#         dict(name='Siemens_SWT_2300kW', up=2500),
#         dict(name='Vestas_V90_3MW', up=50000)
#     ]

#%% Categorise wind turbines function
def categorise_wind_turbines(wind_df,turbine_categories):
    wind_categ = pd.DataFrame(turbine_categories)
    wind_categ['capacity']= 0
    for c in range(len(turbine_categories)):
        if c==0:
            low = 0
        else:
            low = turbine_categories[c-1]['up']
        if c==(len(turbine_categories)-1):
            high = turbine_categories[c]['up']
        else:
            high = turbine_categories[c]['up']
        wind_categ.loc[c,'capacity'] = wind_df.loc[wind_df['capacity'].between(low,high+1,inclusive=False),'capacity'].sum()
    return wind_categ

#%% Categorise wind turbines
# wind_categ = categorise_wind_turbines(onshore_df, turbine_categories)

def capacity_layout(cutout, typ, cap_range=None, until=None):
    """Aggregate selected capacities to the cutouts grid into a capacity layout.

    Parameters
    ----------
        cutout : atlite.cutout
            The cutout for which the capacity layout is contructed.
        typ : str
            Type of energy source, e.g. "Solarstrom" (PV), "Windenergie" (wind).
        cap_range : (optional) list-like
            Two entries, limiting the lower and upper range of capacities (in kW)
            to include. Left-inclusive, right-exclusive.
        until : str
            String representation of a datetime object understood by pandas.to_datetime()
            for limiting to installations existing until this datetime.

    """

    # Load locations of installed capacities and remove incomplete entries
    cols = OrderedDict((('installation_date', 0),
                        ('plz', 2), ('city', 3),
                        ('type', 6),
                        ('capacity', 8), ('level', 9),
                        ('lat', 19), ('lon', 20),
                        ('validation', 22)))
    database = pd.read_csv('eeg_anlagenregister_2015.08.utf8.csv',
                       sep=';', decimal=',', thousands='.',
                       comment='#', header=None,
                       usecols=list(cols.values()),
                       names=list(cols.keys()),
                       # German postal codes can start with '0' so we need to treat them as str
                       dtype={'plz':str},
                       parse_dates=['installation_date'],
                       na_values=('O04WF', 'keine'))

    database = database[(database['validation'] == 'OK') & (database['plz'].notna())]

    # Query postal codes <-> coordinates mapping
    de_nomi = pgeocode.Nominatim('de')
    plz_coords = de_nomi.query_postal_code(database['plz'].unique())
    plz_coords = plz_coords.set_index('postal_code')

    # Fill missing lat / lon using postal codes entries
    database.loc[database['lat'].isna(), 'lat'] = database['plz'].map(plz_coords['latitude'])
    database.loc[database['lon'].isna(), 'lon'] = database['plz'].map(plz_coords['longitude'])

    # Ignore all locations which have not be determined yet
    database = database[database['lat'].notna() & database['lon'].notna()]

    # Select data based on type (i.e. solar/PV, wind, ...)
    data = database[database['type'] == typ].copy()

    # Optional: Select based on installation day
    if until is not None:
        data = data[data['installation_date'] < pd.to_datetime(until)]

    # Optional: Only installations within this caprange (left inclusive, right exclusive)
    if cap_range is not None:
        data = data[(cap_range[0] <= data['capacity']) & (data['capacity'] < cap_range[1])]

    # Determine nearest cells from cutout
    cells = gpd.GeoDataFrame({'geometry': cutout.grid_cells,
                              'lon': cutout.grid_coordinates()[:,0],
                              'lat': cutout.grid_coordinates()[:,1]})

    nearest_cell = cutout.data.sel({'x': data.lon.values,
                                    'y': data.lat.values},
                                   'nearest').coords

    # Map capacities to closest cell coordinate
    data['lon'] = nearest_cell.get('lon').values
    data['lat'] = nearest_cell.get('lat').values

    new_data = data.merge(cells, how='inner')

    # Sum capacities for each grid cell (lat, lon)
    # then: restore lat lon as coumns
    # then: rename and reindex to match cutout coordinates
    new_data = new_data.groupby(['lat','lon']).sum()

    layout = new_data.reset_index().rename(columns={'lat':'y','lon':'x'})\
                        .set_index(['y','x']).capacity\
                        .to_xarray().reindex_like(cutout.data)

    layout = (layout/1e3).fillna(.0).rename('Installed Capacity [MW]')

    return layout