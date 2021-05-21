
from pathlib import Path
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely

def get_countries_regions_ffe(force_recalc=False):
    # Download the region types
    filepath = Path(__file__).parent.parent
    if filepath.joinpath("data_out/zones/zones.csv").is_file() and not force_recalc:
        zones = pd.read_csv(filepath.joinpath("data_out/zones/zones.csv"), index_col=0)
        zones['geometry'] = zones['geometry'].apply(shapely.wkt.loads)
        zones = gpd.GeoDataFrame(zones, geometry="geometry")

    else:
        print("Downloading zone data from FFE ")
        url = "http://opendata.ffe.de:3000/rpc/map_region_type?idregiontype=35"
        zones = gpd.read_file(url)
        zones = zones[['name', 'name_short', 'area_m2', "geometry"]].set_index("name_short")
    
    if filepath.joinpath("data_out/zones/nuts_data.csv").is_file() and not force_recalc:
        nuts_data = pd.read_csv(filepath.joinpath("data_out/zones/nuts_data.csv"), index_col=0)
        nuts_data['geometry'] = nuts_data['geometry'].apply(shapely.wkt.loads)
        nuts_data = gpd.GeoDataFrame(nuts_data, geometry="geometry")

    else:
        print("Downloading nuts data from FFE ")
        url = "http://opendata.ffe.de:3000/rpc/map_region_type?idregiontype=38"
        nuts_data = gpd.read_file(url)
        nuts_data = nuts_data[['id_region', 'name', 'name_short', 'area_m2', "geometry"]]
        nuts_data["country"] = nuts_data.name_short.str[:2]
    return zones, nuts_data


def distance(lat_nodes, lon_nodes, lat_plants, lon_plants):
    """Return vector of distances in km from plant to all nodes in zone."""
    lat_nodes, lon_nodes = lat_nodes.copy(), lon_nodes.copy()
    R = 6373.0

    lat_nodes *= np.pi / 180
    lon_nodes *= np.pi / 180

    lat_plants *= np.pi / 180
    lon_plants *= np.pi / 180

    dlon = lon_nodes - lon_plants
    dlat = lat_nodes - lat_plants

    a = np.sin(dlat / 2)**2 + np.cos(lat_nodes) * np.cos(lat_plants) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R * c

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
    
    zones, nuts_data = get_countries_regions_ffe()
    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    zones.to_csv(wdir.joinpath('data_out/zones/zones.csv'))
    nuts_data.to_csv(wdir.joinpath('data_out/zones/nuts_data.csv'))
    
    # url = "http://opendata.ffe.de:3000/map_region_type?idregiontype=61"
    # zones = gpd.read_file(url)
    # zones = zones[['name', 'name_short', 'area_m2', "geometry"]].set_index("name_short")
    
    
    
    # url = "http://opendata.ffe.de:3000/rpc/region?id_region_type=eq.38"
    # sea_region = gpd.read_file(url)
