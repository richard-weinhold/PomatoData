
from pathlib import Path
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely
import requests 
import json 

import pomato_data

def get_countries_regions_ffe(force_recalc=False):
    # Download the region types

    filepath = Path(pomato_data.__path__[0]).parent  
    if filepath.joinpath("data_out/zones/zones.csv").is_file() and not force_recalc:
        zones = pd.read_csv(filepath.joinpath("data_out/zones/zones.csv"), index_col=0)
        zones['geometry'] = zones['geometry'].apply(shapely.wkt.loads)
        zones = gpd.GeoDataFrame(zones, geometry="geometry")

    else:
        print("Downloading zone data from FFE ")
        url = "http://opendata.ffe.de:3000/rpc/map_region_type?idregiontype=35"
        zones = gpd.read_file(url)
        zones = zones[['name', 'name_short', 'area_m2', "geometry"]].set_index("name_short")
        iso_a3 = pd.read_csv(filepath.joinpath("data_in/regions/country_codes.csv"), index_col=0, usecols=[1,2])
        iso_a3.columns = ["iso_name"]
        iso_a3 = iso_a3.rename({"GB": "UK"})
        zones = pd.merge(zones, iso_a3[["iso_name"]], how="left", left_index=True, right_index=True)
    
    if filepath.joinpath("data_out/zones/nuts_data.csv").is_file() and not force_recalc:
        nuts_data = pd.read_csv(filepath.joinpath("data_out/zones/nuts_data.csv"), index_col=0)
        nuts_data['geometry'] = nuts_data['geometry'].apply(shapely.wkt.loads)
        nuts_data = gpd.GeoDataFrame(nuts_data, geometry="geometry")
    else:
        print("Downloading nuts data from FFE ")
        url = "http://opendata.ffe.de:3000/rpc/map_region_type?idregiontype=38"
        nuts_data = gpd.read_file(url)
        nuts_data = nuts_data[['id_region', 'name', 'name_short', 'area_m2', "geometry"]]
        geometry = []
        for geom in nuts_data['geometry']:
            if geom.is_valid:
                geometry.append(geom)
            else:
                geometry.append(geom.buffer(0))
                
        nuts_data = gpd.GeoDataFrame(nuts_data, geometry=geometry)
        nuts_data["country"] = nuts_data.name_short.str[:2]
        
        coastal = pd.read_csv(filepath.joinpath("data_in/regions/nuts_3_coastal.csv"), index_col=0)
        nuts_data["coastal"] = False
        nuts_data.loc[nuts_data.name_short.isin(coastal.index), "coastal"] = True
    return zones, nuts_data

def get_eez_ffe(force_recalc=False, geometry=True):
    # Download the region types

    filepath = Path(pomato_data.__path__[0]).parent      
    if filepath.joinpath("data_out/zones/eez.csv").is_file() and not force_recalc:
        if geometry:
            eez = pd.read_csv(filepath.joinpath("data_out/zones/eez.csv"), index_col=0)
            eez['geometry'] = eez['geometry'].apply(shapely.wkt.loads)
            eez = gpd.GeoDataFrame(eez, geometry="geometry")
        else:
            eez = pd.read_csv(filepath.joinpath("data_out/zones/eez_wo_geometry.csv"), index_col=0)
        
    else:
        print("Downloading EEZ data from FFE")
        r = requests.get('http://opendata.ffe.de:3000/region?id_region_type=eq.61')
        eez_region = pd.DataFrame.from_dict(json.loads(r.content))
        geometry = []
        for geom in eez_region['geom_4326']:
            tmp_geom =  shapely.wkb.loads(geom, hex=True)
            if tmp_geom.is_valid:
                geometry.append(tmp_geom)
            else:
                geometry.append(tmp_geom.buffer(0))
        eez_region = gpd.GeoDataFrame(eez_region[["id_region", "name"]], geometry=geometry)
        
        eez_region["iso_name"] = eez_region.name.str[:3]
        
        zones, _ = get_countries_regions_ffe()
        zones = zones["iso_name"].reset_index()
        zones.columns = ["zone", "iso_name"]
        
        eez = pd.merge(eez_region, zones, on="iso_name")
        eez["count"] = ""
        for name in eez.name.unique():
            condition = eez.name == name
            eez.loc[condition, "count"] = [str(i + 1) for i in range(sum(condition))]
            
        eez.loc[:, "name"] = eez.name + "_" + eez["count"]
        eez = eez.drop("count", axis=1)
        
    return eez

def distance(lat_nodes, lon_nodes, lat_plants, lon_plants):
    """Return vector of distances in km from plant to all nodes in zone."""
    lat_nodes, lon_nodes = lat_nodes.copy(), lon_nodes.copy()
    R = 6373.0

    lat_nodes *= np.pi / 180
    lon_nodes *= np.pi / 180

    lat_plants *= np.pi / 180
    lon_plants *= np.pi / 180

    delta_lon = lon_nodes - lon_plants
    delta_lat = lat_nodes - lat_plants

    a = np.sin(delta_lat / 2)**2 + np.cos(lat_nodes) * np.cos(lat_plants) * np.sin(delta_lon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R * c

def match_plants_nodes(plants, nodes):
    """Assign nearest node to plants. Subject to penalty depending on voltage level"""
    # plants, nodes = offshore_plants.copy(), offshore_nodes.copy()
    # plants = data.plants.copy()
    # nodes = data.nodes.copy()
    # p ="H978"
    # plants.loc[p, "node"] = None
    
    grid_node = []
    condition = plants.node.isna()
    for p in plants[condition].index:
        # 1
        nodes_in_area = nodes[nodes.zone == plants.zone[p]].copy()
        if "n" + plants.zone[p] in nodes_in_area.index:
            grid_node.append("n" + plants.zone[p])
        else:
            nodes_in_area["distance"] = distance(
                nodes_in_area.lat.values, nodes_in_area.lon.values, 
                plants.loc[p, "lat"], plants.loc[p, "lon"])
            
            
            if plants.loc[p, "g_max"] > 1000:
                weighting = {500: 1, 380: 1.2, 220: 100, 132: 1000}
                nodes_in_area["distance_penalty"] = nodes_in_area["voltage"].map(weighting)
            elif plants.loc[p, "g_max"] > 500:
                weighting = {500: 1, 380: 1.2, 220: 2, 132: 5}
                nodes_in_area["distance_penalty"] = nodes_in_area["voltage"].map(weighting)
            else:
                weighting = {500: 1, 380: 1.2, 220: 1.3, 132: 1.4}
                nodes_in_area["distance_penalty"] = nodes_in_area["voltage"].map(weighting)
                
            nodes_in_area.loc[:, "distance"] *= nodes_in_area["distance_penalty"]
            
            grid_node.append(nodes_in_area["distance"].idxmin())
            
            if nodes_in_area["distance"].min() > 100 and len(nodes_in_area) > 10:
                zone = nodes.loc[grid_node[-1], "zone"]
                print(f"Plant {p} is more than 100 km away from node {grid_node[-1]} in zone {zone}")
            
    plants.loc[condition.values, "node"] = grid_node  
    
    return plants 

def load_data_structure(wdir):
    """Read in datastructure and init attributes with existing data."""
    file = wdir.joinpath("data_in/data_structure.xlsx")
    xls = pd.ExcelFile(file)
    structure = xls.parse("raw")
    columns = [c for c in structure.columns if "Unnamed:" not in c]
    
    data_structure = {}

    for c in columns:
        att = "attributes"
        col_pos = structure.columns.get_loc(c)
        cols = list(structure.columns[col_pos:col_pos + 2])
        tmp = structure.loc[1:, cols].copy()
        data_structure[c] = {"attributes": {}, "optional attributes": {}}
        for (t, v) in zip(tmp[cols[0]].astype(str), tmp[cols[1]]):
            if not t == "nan":
                if t == "optional attributes":
                    att = "optional attributes"
                else:
                    data_structure[c][att][t] = v
    return data_structure

def add_timesteps(df):
    
    if not ("utc_timestamp" in df.columns or df.index.name == "utc_timestamp"):
        raise ValueError("utc_timestamp not in DataFrame, cannot add timsteps")
    
    if df.index.name == "utc_timestamp":
        df = df.reset_index()
    timestamps = df.utc_timestamp.sort_values().unique()
    timestamps_t_strings = ['t' + "{0:0>4}".format(x) for x in range(1, len(timestamps) + 1)] 
    timestamp_dict = {k: v for k, v in zip(timestamps, timestamps_t_strings)}
    df["timestep"] = df["utc_timestamp"].replace(timestamp_dict)
    if (len(df.index) == len(timestamps)):
        df.set_index("timestep", inplace=True)
    return df

# %%

if __name__ == "__main__":
    wdir = Path(pomato_data.__path__[0]).parent  

    zones, nuts_data = get_countries_regions_ffe(force_recalc=True)
    # zones.to_csv(wdir.joinpath('data_out/zones/zones.csv'))
    # nuts_data.to_csv(wdir.joinpath('data_out/zones/nuts_data.csv'))
    
    # eez_region = get_eez_ffe(force_recalc=True)
    # eez = get_eez_ffe()
    # eez.drop("geometry", axis=1).to_csv(wdir.joinpath('data_out/zones/eez_wo_geometry.csv'))
    # eez.to_csv(wdir.joinpath('data_out/zones/eez.csv'))
    
    # url = "http://opendata.ffe.de:3000/rpc/region?id_region_type=eq.38"
    # sea_region = gpd.read_file(url)
