
import os
import requests
import pandas as pd
from pathlib import Path
import shapely
import geopandas as gpd

os.chdir(r'C:\Users\riw\Documents\repositories\pomato_data')
from pomato_data.auxiliary import get_countries_regions_ffe

def anymod_installed_capacities(base_path, year):
    # base_path = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    anymod_result_path = base_path.joinpath("data_in/anymod_results/results_summary_8days_grid_202105061657.csv")

    raw = pd.read_csv(anymod_result_path)
    raw.loc[raw.technology.str.contains("offshore"), "group"] +=  "offshore" 
    raw.loc[raw.technology.str.contains("onshore"), "group"] += "onshore" 
    raw.loc[:, "group"] = raw.loc[:, "group"].str.rstrip(" ")
    raw = raw.drop(["id", "Unnamed: 10"], axis=1)
    
    installed_capacity  = (raw[(raw.variable.isin(["capaConv", "capaStSize","capaStIn"])) & \
                              # (raw.country == "DE") & \
                              (raw.timestep_superordinate_dispatch.str.contains(str(year)))
                              ].groupby(["country", "group", "variable"]).sum()
                           ).reset_index().pivot(index=["country", "group"], 
                                                 columns="variable", 
                                                 values="value")
    return installed_capacity


def calculate_capacities_from_pontentials(base_path, year=2030):
    # c = "CH"
    # base_path = wdir

    wind_potentials = pd.read_csv(base_path.joinpath('data_out/res_potential/wind_potential.csv'), index_col=0).set_index("name_short")
    pv_potentials = pd.read_csv(base_path.joinpath('data_out/res_potential/pv_potential.csv'), index_col=0).set_index("name_short")
    installed_capacities = anymod_installed_capacities(base_path, year)
    # Capacity is distributed based on the potential in MWh
    wind_capacities = wind_potentials.copy()
    pv_capacities = pv_potentials.copy()
    countries = [c for c in wind_capacities.country.unique() if c in installed_capacities.index]

    wind_capacities["relative_potential"] = 1 
    wind_capacities["capacity"] = 0 
    pv_capacities["relative_potential"] = 1 
    pv_capacities["capacity"] = 0 
    for c in countries:
        wind_capacities.loc[wind_capacities.country == c, "relative_potential"] = wind_potentials.loc[wind_potentials.country == c, "value"]/wind_potentials.loc[wind_potentials.country == c, "value"].sum()
        wind_capacities.loc[wind_capacities.country == c, "capacity"] = wind_capacities.loc[wind_capacities.country == c, "relative_potential"]*installed_capacities.loc[(c, "wind onshore"), "capaConv"]*1000
        pv_capacities.loc[pv_capacities.country == c, "relative_potential"] = pv_potentials.loc[pv_potentials.country == c, "value"]/pv_potentials.loc[pv_potentials.country == c, "value"].sum()
        pv_capacities.loc[pv_capacities.country == c, "capacity"] = pv_capacities.loc[pv_capacities.country == c, "relative_potential"]*installed_capacities.loc[(c, "solar"), "capaConv"]*1000
    return wind_capacities, pv_capacities
    


def regionalize_wind_solar_capacities(wdir, nodes, zones):
    
    capacity_wind = pd.read_csv(wdir.joinpath('data_out/res_capacity/wind_capacity.csv'), 
                                index_col=0).rename(columns={"capacity": "wind onshore"})
    capacity_pv = pd.read_csv(wdir.joinpath('data_out/res_capacity/pv_capacity.csv'), 
                              index_col=0).rename(columns={"capacity": "solar"})
    
    capacity_nuts = pd.merge(capacity_pv["solar"], capacity_wind["wind onshore"], 
                             right_index=True, left_index=True)
    
    plants = pd.DataFrame()
    for zone in zones:
        plants = pd.concat([plants, regionalize_capacities_country(nodes, zone, capacity_nuts)])
    return plants
        
def regionalize_capacities_country(nodes, zone, capacity_nuts):
    # nodes = self.nodes.copy()
    # # nodes[nodes.zone == zone]
    # zone = "NO"    
    print("Regionalizing wind/pv capacities for", zone)
    
    nodes = nodes[nodes["info"] != "joint"].copy()
    country_data, nuts_data = get_countries_regions_ffe()    
    # country_data['geometry'] = country_data['geometry'].apply(shapely.wkt.loads)
    # nuts_data['geometry'] = nuts_data['geometry'].apply(shapely.wkt.loads)
    
    nuts_data = gpd.GeoDataFrame(nuts_data[nuts_data.country == zone]).set_crs("EPSG:4326")

    geometry = [shapely.geometry.Point(xy) for xy in zip(nodes.lon, nodes.lat)]
    cond = [n.within(country_data.loc[zone, "geometry"]) for n in geometry]
    nodes = gpd.GeoDataFrame(nodes, crs="EPSG:4326", geometry=geometry).loc[cond]
    
    capacity_nodes = nodes[['voltage', 'name', 'zone', "info", "lat", "lon"]].copy()
    nuts_to_nodes = gpd.sjoin(nodes, nuts_data, how='left', op='within')

    nuts_to_nodes["weighting"] = 1
    nuts_multiple_nodes = nuts_to_nodes[["name_short", "weighting"]].groupby('name_short').count()
    nuts_multiple_nodes["sum_weight"] = nuts_to_nodes[["name_short", "weighting"]].groupby('name_short').sum()
    
    capacity_nodes["wind onshore"] = (capacity_nuts.loc[nuts_to_nodes.name_short, "wind onshore"]/nuts_multiple_nodes.loc[nuts_to_nodes.name_short, "sum_weight"]).values
    capacity_nodes["solar"] = (capacity_nuts.loc[nuts_to_nodes.name_short, "solar"]/nuts_multiple_nodes.loc[nuts_to_nodes.name_short, "sum_weight"]).values
   
    node_in_nuts = gpd.sjoin(nuts_data, nodes, how='left', op='contains')
    nuts_no_node = node_in_nuts[node_in_nuts["index_right"].isna()]
    
    nuts_centroids = nuts_data.centroid # .to_crs("epsg:4326")
    nuts_no_node = nuts_centroids.loc[nuts_no_node.index] #.to_crs('epsg:2953')

    # Find the closest node to the region's centroid and map them
    dist = []
    for j in range(len(nodes)):
        dist.append(nuts_no_node.distance(nodes.geometry.iloc[j]))
    dist = pd.DataFrame(dist)
        
    nuts_no_node_map = pd.DataFrame(dist.idxmin(), columns=['node_index'])
    nuts_no_node_map["node"] = nodes.iloc[nuts_no_node_map.node_index].index
    nuts_no_node_map = nuts_no_node_map.join(pd.DataFrame(nuts_data['name_short']))
    nuts_no_node_map = nuts_no_node_map.reset_index()
    nuts_no_node_map.rename(columns={'index': 'nuts_index'}, inplace=True)
    
    # Add those regions to the nodal load patterns
    for i in nuts_no_node_map.index:
        capacity_nodes.loc[nuts_no_node_map['node'].loc[i], "wind onshore"] += capacity_nuts.loc[nuts_no_node_map['name_short'].loc[i], "wind onshore"]
        capacity_nodes.loc[nuts_no_node_map['node'].loc[i], "solar"] += capacity_nuts.loc[nuts_no_node_map['name_short'].loc[i], "solar"]

    wind_plants = capacity_nodes[["zone", "wind onshore", 'name', "lat", "lon"]].reset_index()
    wind_plants.columns = ["node", "zone", "g_max", 'name', "lat", "lon"]
    wind_plants["index"] = wind_plants.node.astype(str) + "_wind"
    wind_plants = wind_plants.set_index("index")
    wind_plants["eta"] = 1
    wind_plants["h_max"] = 0
    wind_plants["fuel"] = "wind"
    wind_plants["commissioned"] = 1900
    wind_plants["status"] = "online"
    wind_plants[["plant_type", "technology"]] = "wind onshore"
    wind_plants[["chp", "heatarea", "city", "postcode", "street"]] = None
    
    pv_plants = capacity_nodes[["zone", "solar", 'name', "lat", "lon"]].reset_index()
    pv_plants.columns = ["node", "zone", "g_max", 'name', "lat", "lon"]
    pv_plants["index"] = pv_plants.node.astype(str) + "_solar"
    pv_plants = pv_plants.set_index("index")
    pv_plants["eta"] = 1
    pv_plants["h_max"] = 0
    pv_plants["fuel"] = "sun"
    pv_plants["commissioned"] = 1900
    pv_plants["status"] = "online"
    pv_plants[["plant_type", "technology"]] = "solar"
    pv_plants[["chp", "heatarea", "city", "postcode", "street"]] = None
    
    return pd.concat([wind_plants, pv_plants])


# %%
if __name__ == "__main__": 
    
    from shapely import wkt

    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    wind_capacities, pv_capacities = calculate_capacities_from_pontentials(wdir)
    wind_capacities.to_csv(wdir.joinpath('data_out/res_capacity/wind_capacity.csv'))
    pv_capacities.to_csv(wdir.joinpath('data_out/res_capacity/pv_capacity.csv'))

    # country_data, nuts_data = get_countries_regions_ffe()    
    # df = pd.merge(wind_capacities, nuts_data, left_index=True, right_on="name_short")
    wind_capacities['geometry'] = wind_capacities['geometry'].apply(wkt.loads)
    wind_capacities.loc[wind_capacities.capacity > 2000, "capacity"] = 2000
    gpd.GeoDataFrame(wind_capacities, geometry="geometry").plot(column="capacity", legend=True)  
    pv_capacities['geometry'] = pv_capacities['geometry'].apply(wkt.loads)
    gpd.GeoDataFrame(pv_capacities, geometry="geometry").plot(column="capacity", legend=True)  

    installed_capacities = anymod_installed_capacities(wdir, 2030)

    
    