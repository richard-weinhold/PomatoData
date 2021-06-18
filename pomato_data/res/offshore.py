
import os
import requests
import pandas as pd
from pathlib import Path
import shapely
import pyproj
import geopandas as gpd

from shapely.geometry import Point, LineString

from pomato_data.auxiliary import get_countries_regions_ffe, match_plants_nodes, get_eez_ffe
from pomato_data.res import anymod_installed_capacities


def process_offshore_windhubs(wdir, nodes, weather_year, capacity_year):
    # weather_year, capacity_year = 2019, 2020
    
    offshore_plants, offshore_nodes = create_offshore_hubs(wdir, weather_year, capacity_year)

    offshore_connections = pd.read_csv(wdir.joinpath("data_in/nodes/offshore_connections.csv"))
    offshore_connections = offshore_connections[offshore_connections.node.isin(nodes.index)]
    zones = nodes.zone.unique()
    
    offshore_plants = offshore_plants[offshore_plants.zone.isin(zones)]
    offshore_nodes = offshore_nodes[offshore_nodes.index.isin(offshore_plants["node"])]
    
    geod = pyproj.Geod(ellps="WGS84")
    add_dclines_data = []
    add_dclines_index = []
    for node_j in offshore_nodes.index:
        onshore_nodes = offshore_connections.loc[offshore_connections.offshore_hub == node_j, "node"]
        if onshore_nodes.empty:
            zone = offshore_nodes.loc[node_j, "zone"]
            node_i = "n" + zone
            geometry = LineString([Point(nodes.loc[node_i, ["lon", "lat"]]), 
                                   Point(offshore_nodes.loc[node_j, ["lon", "lat"]])])
            length = geod.geometry_length(geometry)
            add_dclines_index.append("dc" + node_j)
            add_dclines_data.append([node_i, node_j, zone, node_j, 500, length, False,
                                geometry.wkt, "dc", 1e5, "online", 1900])
            
        elif all(n in nodes.index for n in onshore_nodes):
            for node_i in onshore_nodes:
                geometry = LineString([Point(nodes.loc[node_i, ["lon", "lat"]]), 
                                       Point(offshore_nodes.loc[node_j, ["lon", "lat"]])])
                length = geod.geometry_length(geometry)
                zone = nodes.loc[node_i, "zone"]
                add_dclines_index.append("dc" + node_j + "_" + node_i)
                add_dclines_data.append([node_i, node_j, zone, node_j, 500, length, False,
                                         geometry.wkt, "dc", 1e5, "online", 1900])
        
    cols = ['node_i', 'node_j', 'name_i', 'name_j', 'voltage', 'length', 'under_construction', 
            'geometry', 'technology', 'capacity', 'status', 'commissioned']
    dclines = pd.DataFrame(index=add_dclines_index, columns=cols,
                           data=add_dclines_data)
    
    eez = get_eez_ffe(geometry=False).set_index("name")
    
    availability_raw = pd.read_csv(wdir.joinpath(f'data_out/res_availability/offshore_availability_{weather_year}.csv'))
    availability_raw = availability_raw.pivot(index="utc_timestamp", columns="id_region", values="value")
    availability_raw.index = pd.to_datetime(availability_raw.index).astype('datetime64[ns]')
    
    availability_raw = availability_raw.sort_index()
    
    availability = pd.DataFrame(index=availability_raw.index)
    for p in offshore_plants.index:
        availability[p] = availability_raw[eez.loc[p.lstrip("n"), "id_region"]]
    
    return offshore_plants, offshore_nodes, dclines, availability

def create_offshore_hubs(wdir, weather_year, capacity_year):
    # capacity_year = 2019
    eez = get_eez_ffe()
    # eez = eez.set_index("name")
    # eez.loc[["BEL_north_sea_1"], :].plot()
    
    country_data, nuts_data = get_countries_regions_ffe()    
    
    offshore_nodes_index, offshore_nodes_data = [], [] 
    for z in eez.index:
        name = eez.loc[z, "name"]
        zone = eez.loc[z, "zone"]
        idx = "n" + name
        lon, lat = eez.loc[z, "geometry"].centroid.coords[0]
        offshore_nodes_data.append([idx, 500, name, lat, lon, zone, "", False, True])
        offshore_nodes_index.append(idx)
    cols = ['substation', 'voltage', 'name', 'lat', 'lon', 'zone', 'info', 'demand', 'slack']
    offshore_nodes = pd.DataFrame(index=offshore_nodes_index, columns=cols,
                                  data=offshore_nodes_data)
    
    corrd_correction = {
        "nDEU_baltic_sea_1": (54.753567, 13.974048),
        "nDNK_baltic_sea_1": (56.574315, 11.634032),
        "nDNK_north_sea_1": (56.425272, 7.825166),
        "nDNK_north_sea_1": (56.425272, 7.825166),
        "nFRA_north_sea_1": (50.095828, -0.238387),
        "nSWE_baltic_sea_1": (62.427037, 19.507758),
        "nESP_atlantic_1": (44.370247, -6.348793),
        "nESP_atlantic_2": (36.415398, -7.185952),
        }
    
    for n in corrd_correction:
        offshore_nodes.loc[n, ["lat", "lon"]] = corrd_correction[n]
    
    installed_capacities = anymod_installed_capacities(wdir, capacity_year)
    installed_capacities.xs("wind offshore", level=1)
    installed_capacities.xs("stock", level=1)
    
    offshore_plants = offshore_nodes[['name', 'lat', 'lon', 'zone']].copy()
    eez.loc[:, "name"] = "n" +  eez.name
    
    availability_factor = pd.read_csv(wdir.joinpath(f'data_out/res_availability/offshore_availability_{weather_year}.csv'), index_col=0).groupby("id_region").sum()
    
    offshore_plants["g_max"] = 0
    for zone in offshore_plants.zone.unique():
        if zone in installed_capacities.xs("wind offshore", level=1).index:
            # zone = "DE"
            capacity = installed_capacities.loc[(zone, "wind offshore"), "capaConv"]*1000
            
            tmp = eez.loc[eez.zone == zone, ["id_region", "name"]].set_index("id_region")
            
            tmp["g_max"] = availability_factor.loc[tmp.index] / availability_factor.loc[tmp.index].sum() * capacity
            offshore_plants.loc[tmp.name, "g_max"] = tmp.g_max.values
            
        else:
            offshore_plants = offshore_plants[offshore_plants.zone != zone]
            
    offshore_plants = offshore_plants[offshore_plants.g_max > 0]    
    offshore_plants["technology"] = "wind offshore"
    offshore_plants["plant_type"] = "wind offshore"
    offshore_plants["fuel"] = "wind"
    offshore_plants["eta"] = 1
    offshore_plants[["mc_el", "mc_heat"]] = 5, 0
    offshore_plants["commissioned"] = 1900
    offshore_plants["status"] = "online"
    offshore_plants[["chp", "heatarea", "city", "postcode", "street", "company", "storage_capacity"]] = None
    offshore_plants["node"] = offshore_plants.index
    offshore_plants[["d_max", "h_max"]] = 0
    
    cols = ['g_max', 'chp', 'h_max', 'city', 'commissioned', 'company', 'zone',
            'eta', 'fuel', 'lat', 'lon', 'name', 'postcode', 'status', 'street',
            'technology', 'plant_type', 'heatarea', 'storage_capacity', 'node',
            'mc_el', 'mc_heat', 'd_max']
    
    return offshore_plants[cols], offshore_nodes
    
# %%
if __name__ == "__main__": 
    
    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    offshore_plants, offshore_nodes = create_offshore_hubs(wdir, 2019, 2020)
    # offshore_plants.to_csv(wdir.joinpath("data_out/res_capacity/offshore.csv"))
    # offshore_nodes.to_csv(wdir.joinpath("data_out/res_capacity/offshore_nodes.csv"))
    
    
    
    