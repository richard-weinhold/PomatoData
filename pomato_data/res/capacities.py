
import pandas as pd
from pathlib import Path
import shapely
import geopandas as gpd

from pomato_data.auxiliary import get_countries_regions_ffe, match_plants_nodes, get_eez_ffe

def anymod_installed_capacities(wdir, year=2020):
    # base_path = Path(pomato_data.__path__[0]).parent 
    anymod_result_path = wdir.joinpath("data_in/anymod_results/results_summary_8days_grid_202105061657.csv")

    raw = pd.read_csv(anymod_result_path)
    raw.loc[raw.technology.str.contains("offshore"), "group"] +=  "offshore" 
    raw.loc[raw.technology.str.contains("onshore"), "group"] += "onshore" 
    raw.loc[:, "group"] = raw.loc[:, "group"].str.rstrip(" ")
    raw = raw.drop(["id", "Unnamed: 10"], axis=1)
    
    condition_var = raw.variable.isin(["capaConv", "capaStSize","capaStIn"])
    condition_year = raw.timestep_superordinate_dispatch.str.contains(str(year))
    
    if sum(condition_year) == 0:
        raise ValueError("capacity year not in AnyMod results")
    else:
        installed_capacity = raw[condition_var & condition_year].groupby(["country", "group", "variable"]).sum()
        installed_capacity  = installed_capacity.reset_index().pivot(
            index=["country", "group"],columns="variable", values="value")
                                     
    installed_capacity = installed_capacity[installed_capacity.index.isin(["solar", "wind onshore", "wind offshore"], level=1)]
    return installed_capacity

def input_installed_capacities(wdir, filepath):
    filepath = Path(r"C:\Users\riw\Documents\repositories\pomato_data\data_in\res\manual_capacities_2022.csv")
    capacities = pd.read_csv(filepath, index_col=[0,1])
    return capacities 


def calculate_capacities_from_potentials(wdir, settings):

    wind_potentials = pd.read_csv(wdir.joinpath('data_out/res_potential/wind_potential.csv'), index_col=0).set_index("name_short")
    pv_potentials = pd.read_csv(wdir.joinpath('data_out/res_potential/pv_potential.csv'), index_col=0).set_index("name_short")
    
    if settings["capacity_source"] == "anymod":
        installed_capacities = anymod_installed_capacities(wdir, settings["capacity_year"])
    else:
        installed_capacities = input_installed_capacities(wdir, settings["capacity_file"])

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
        if c in installed_capacities.xs("wind onshore", level=1).index:
            wind_capacities.loc[wind_capacities.country == c, "capacity"] = wind_capacities.loc[wind_capacities.country == c, "relative_potential"]*installed_capacities.loc[(c, "wind onshore"), "capaConv"]*1000
        else:
            wind_capacities.loc[wind_capacities.country == c, "capacity"] = 0
        
        pv_capacities.loc[pv_capacities.country == c, "relative_potential"] = pv_potentials.loc[pv_potentials.country == c, "value"]/pv_potentials.loc[pv_potentials.country == c, "value"].sum()
        if c in installed_capacities.xs("solar", level=1).index:
            pv_capacities.loc[pv_capacities.country == c, "capacity"] = pv_capacities.loc[pv_capacities.country == c, "relative_potential"]*installed_capacities.loc[(c, "solar"), "capaConv"]*1000
        else:
            pv_capacities.loc[pv_capacities.country == c, "capacity"] = 0

    return wind_capacities, pv_capacities
    
def regionalize_res_capacities(wdir, nodes, zones, technology, settings):
        
    capacity_wind, capacity_pv = calculate_capacities_from_potentials(wdir, settings)
    
    capacity_wind = capacity_wind.rename(columns={"capacity": "wind/wind onshore"})
    capacity_pv = capacity_pv.rename(columns={"capacity": "sun/solar"})
    
    other_res = pd.read_csv(wdir.joinpath('data_out/res_capacity/other_res.csv'), 
                              index_col=0)
    
    capacity_nuts = pd.merge(capacity_pv["sun/solar"], capacity_wind["wind/wind onshore"], 
                             right_index=True, left_index=True)
    capacity_nuts = pd.merge(capacity_nuts, other_res, right_index=True, left_index=True, how="left").fillna(0)
    
    plants = pd.DataFrame()
    for zone in zones:
        plants = pd.concat([plants, regionalize_capacities_country(nodes, zone, capacity_nuts, technology)])
    return plants
        
def regionalize_capacities_country(nodes, zone, capacity_nuts, technology):
    # nodes = data.nodes.copy()
    # technology = data.technology.copy()
    # # # nodes[nodes.zone == zone]
    # zone = "DE"    
    print("Regionalizing wind/pv capacities for", zone)
    
    nodes = nodes[nodes["info"] != "joint"].copy()
    country_data, nuts_data = get_countries_regions_ffe()    
    # country_data['geometry'] = country_data['geometry'].apply(shapely.wkt.loads)
    # nuts_data['geometry'] = nuts_data['geometry'].apply(shapely.wkt.loads)
    
    nuts_data = gpd.GeoDataFrame(nuts_data[nuts_data.country == zone]).set_crs("EPSG:4326")

    geometry = [shapely.geometry.Point(xy) for xy in zip(nodes.lon, nodes.lat)]
    condition = [n.within(country_data.loc[zone, "geometry"]) for n in geometry]
    nodes = gpd.GeoDataFrame(nodes, crs="EPSG:4326", geometry=geometry).loc[condition]
    
    capacity_nodes = nodes[['voltage', 'name', 'zone', "info", "lat", "lon"]].copy()
    nuts_to_nodes = gpd.sjoin(nodes, nuts_data, how='left', op='within')

    nuts_to_nodes["weighting"] = 1
    nuts_multiple_nodes = nuts_to_nodes[["name_short", "weighting"]].groupby('name_short').count()
    nuts_multiple_nodes["sum_weight"] = nuts_to_nodes[["name_short", "weighting"]].groupby('name_short').sum()

    for tech in capacity_nuts.columns:
        capacity_nodes[tech] = (capacity_nuts.loc[nuts_to_nodes.name_short, tech]/nuts_multiple_nodes.loc[nuts_to_nodes.name_short, "sum_weight"]).values
   
    node_in_nuts = gpd.sjoin(nuts_data, nodes, how='left', op='contains')
    nuts_no_node = node_in_nuts[node_in_nuts["index_right"].isna()]
    
    nuts_centroids = nuts_data.centroid # .to_crs("epsg:4326")
    # nuts_centroids = nuts_data.geometry.to_crs("epsg:3395").centroid.to_crs('epsg:4326')
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
        for tech in capacity_nuts.columns:
            capacity_nodes.loc[nuts_no_node_map['node'].loc[i], tech] += capacity_nuts.loc[nuts_no_node_map['name_short'].loc[i], tech]


    plants = pd.DataFrame()
    technology = technology.set_index(["fuel", "technology"])
    for tech in capacity_nuts.columns:
        tmp_plants = capacity_nodes[["zone", 'name', "lat", "lon", tech]].reset_index().rename(columns={tech: "g_max", "index": "node"})
        tmp_plants["index"] = tmp_plants.node.astype(str) + "_" + tech
        tmp_plants = tmp_plants.set_index("index")
        tmp_plants[["fuel", "technology"]] = tech.split("/")
        tmp_plants["eta"] = technology.loc[tuple(tech.split("/")), "eta"]
        tmp_plants["plant_type"] = technology.loc[tuple(tech.split("/")), "plant_type"]
        tmp_plants["commissioned"] = 1900
        tmp_plants["status"] = "online"
        tmp_plants[["chp", "heatarea", "city", "postcode", "street"]] = None
        tmp_plants[tmp_plants.g_max > 0]
        plants = pd.concat([plants, tmp_plants[tmp_plants.g_max > 0]])   
    
    return plants

def existing_offshore_wind_capacities(wdir, nodes, weather_year):
    
    plants = pd.read_csv(wdir.joinpath('data_in/res/renewable_power_plants_EU.csv'))
    condition = (plants.technology == "Offshore")&(~plants.lat.isna()&(plants.country.isin(nodes.zone)))
    offshore_plants = plants.loc[condition]
    cols = ["country", "lat", "lon", "electrical_capacity"]
    
    offshore_plants = offshore_plants[cols]
    new_cols = ["zone", "lat", "lon", "g_max"]
    offshore_plants.columns =  new_cols
    
    offshore_plants["technology"] = "wind offshore"
    offshore_plants["plant_type"] = "wind offshore"
    offshore_plants["fuel"] = "wind"
    offshore_plants["eta"] = 1
    offshore_plants[["mc_el", "mc_heat"]] = 0,0
    offshore_plants["commissioned"] = 1900
    offshore_plants["status"] = "online"
    offshore_plants[["chp", "heatarea", "city", "postcode", "street"]] = None
    offshore_plants["node"] = None

    eez = get_eez_ffe()
    eez = eez.set_index("id_region").rename(columns={"name": "name_eez"})

    offshore_nodes = nodes.copy()
    geometry = [shapely.geometry.Point(xy) for xy in zip(offshore_nodes.lon, offshore_nodes.lat)]
    offshore_nodes = gpd.GeoDataFrame(offshore_nodes, crs="EPSG:4326", geometry=geometry)
    
    offshore_nodes = gpd.sjoin(offshore_nodes, eez[["name_eez", "geometry"]].set_crs(crs="EPSG:4326"), how='left', op='within')
    
    offshore_nodes = offshore_nodes[(~offshore_nodes.name_eez.isna())|(offshore_nodes.voltage >= 500)]
    offshore_plants = match_plants_nodes(offshore_plants, offshore_nodes)
    offshore_plants = offshore_plants.reset_index()
    offshore_plants["index"] = offshore_plants.node.astype(str) + "_offshore" + "_" + offshore_plants.index.astype(str)
    offshore_plants = offshore_plants.set_index("index")
    
    geometry = [shapely.geometry.Point(xy) for xy in zip(offshore_plants.lon, offshore_plants.lat)]
    offshore_plants = gpd.GeoDataFrame(offshore_plants, crs="EPSG:4326", geometry=geometry)
    eez_centroids = eez.centroid # .to_crs("epsg:4326")
    dist = []
    for j in range(len(offshore_plants)):
        dist.append(eez_centroids.distance(offshore_plants.geometry.iloc[j]))
    dist = pd.DataFrame(dist)
        
    nuts_no_node_map = pd.DataFrame(dist.idxmin(axis=1), columns=['eez_index'])
    offshore_plants["eez_id"] = nuts_no_node_map.eez_index.values
    
    availability_raw = pd.read_csv(wdir.joinpath(f'data_out/res_availability/offshore_availability_{weather_year}.csv'))
    availability_raw = availability_raw.pivot(index="utc_timestamp", columns="id_region", values="value")
    availability_raw.index = pd.to_datetime(availability_raw.index).astype('datetime64[ns]')
    availability_raw = availability_raw.sort_index()
    
    availability = pd.DataFrame(index=availability_raw.index)
    for p in offshore_plants.index:
        availability[p] = availability_raw[offshore_plants.loc[p, "eez_id"]]
    
    return availability, offshore_plants
    
def other_res_capacities(wdir):
    
    plants_raw = pd.read_csv(wdir.joinpath('data_in/res/renewable_power_plants_EU.csv'))
    condition = (plants_raw.energy_source_level_2 != "Wind")&(plants_raw.energy_source_level_2 != "Solar")
    
    plants = plants_raw.loc[condition]
    condition = plants.nuts_3_region.isna()
    
    condition = plants.lat.isna()
    remove_capacity = plants.loc[condition, "electrical_capacity"].sum()
    print(f"Remocing {remove_capacity.round()} MW of capacity because no NUTS3 information is available")
    
    plants = plants[~condition]
    plants = plants.rename(columns={"energy_source_level_3": "fuel"})
    plants.loc[plants.fuel.isna(), "fuel"] = plants.loc[plants.fuel.isna(), "energy_source_level_2"]
    
    fuel_dict = {"Hydro": "hydro", 
                 "Geothermal": "sun", 
                 'Biomass and biogas': "biomass", "Biomass and Biogas": "biomass",
                 "Other or unspecified": "biomass", "Sewage and landfill gas": "biomass",
                 "Other bioenergy and renewable waste": "biomass"}
    
    tech_dict = {'Run-of-river': "ror", "Pumped storage": "psp", "Combustion engine": "generator",
                 "sewage and landfill gas": "ccgt", "Sewage and landfill gas": "ccgt",
                 "Biomass and biogas": "ccgt",
                 "Geothermal": "geothermal", "Other or unspecified": "biomass"}
    
    plants.technology = plants.technology.replace(tech_dict, regex=False)
    plants.fuel = plants.fuel.replace(fuel_dict, regex=False)
    
    plants.loc[(plants.fuel == "Marine"), ["technology", "fuel"]] = "tidal", "hydro"
    plants.loc[(plants.technology == "Other or unspecified technology")&(plants.fuel == "biomass"), "technology"] = "ccgt"
    plants.loc[(plants.technology.isna())&(plants.fuel == "biomass"), "technology"] = "ccgt"
    
    plants.loc[(plants.technology == "Other or unspecified technology")&(plants.fuel == "hydro"), "technology"] = "hydro"
    plants.loc[(plants.technology == "Other or unspecified technology")&(plants.energy_source_level_2 == "Hydro"), ["technology", "fuel"]] = "hydro", "hydro"
    plants.loc[(plants.technology.isna())&(plants.fuel == "hydro"), "technology"] = "hydro"
    
    plants.loc[(plants.energy_source_level_2 == "Geothermal"), "technology"] = "geothermal"

    cols = ["fuel", "technology", "nuts_3_region", "electrical_capacity"]
    
    plants = plants[cols].groupby(cols[:-1]).sum().reset_index()
    plants.columns = ["fuel", "technology", "nuts_id", "capacity"]
    
    return plants

# %%
if __name__ == "__main__": 
    import pomato_data
    import matplotlib.pyplot as plt        

    wdir = Path(pomato_data.__path__[0]).parent 
    wind_capacities, pv_capacities = calculate_capacities_from_potentials(wdir, 2030)
    
    wind_capacities.loc[wind_capacities.country == "DE", "capacity"].sum()
    pv_capacities.loc[pv_capacities.country == "DE", "capacity"].sum()
    
    wind_capacities = wind_capacities.loc[wind_capacities.country == "DE", :]
    pv_capacities = pv_capacities.loc[pv_capacities.country == "DE", :]
    # wind_capacities.to_csv(wdir.joinpath('data_out/res_capacity/wind_capacity.csv'))
    # pv_capacities.to_csv(wdir.joinpath('data_out/res_capacity/pv_capacity.csv'))
    
    # other_res_capacities = other_res(wdir)
    # other_res_capacities.to_csv(wdir.joinpath('data_out/res_capacity/other_res.csv'))

    country_data, nuts_data = get_countries_regions_ffe()    
    df = pd.merge(wind_capacities, nuts_data, left_index=True, right_on="name_short")
    
    filepath = Path(r"C:\Users\riw\tubCloud\Uni\Projekte\Modezeen\Arbeitstreffen 2\plots")

    wind_capacities['geometry'] = wind_capacities['geometry'].apply(shapely.wkt.loads)
    # wind_capacities.loc[wind_capacities.capacity > 2000, "capacity"] = 1800
    # wind_capacities.loc[wind_capacities.capacity > wind_capacities.capacity.quantile(0.95), "capacity"] = wind_capacities.capacity.quantile(0.95)
    wind_capacities["capacity_per_km2"] = wind_capacities.capacity / (wind_capacities.area_m2 / 1e6)

    fig, ax = plt.subplots()
    gpd.GeoDataFrame(wind_capacities, geometry="geometry").plot(
        column="capacity", legend=True, cmap="summer", ax=ax,
        legend_kwds={"label": "Installed Capacity in MW"})  
    fig.savefig(filepath.joinpath("wind_de.png"))

    fig, ax = plt.subplots()
    gpd.GeoDataFrame(wind_capacities, geometry="geometry").plot(
        column="capacity_per_km2", legend=True, cmap="summer", ax=ax,
        legend_kwds={"label": "Installed Capacity per $km^2$ in MW"}) 
    fig.savefig(filepath.joinpath("wind_de_per_km2.png"))

    pv_capacities.loc[:, 'geometry'] = pv_capacities['geometry'].apply(shapely.wkt.loads)
    pv_capacities["capacity_per_km2"] = pv_capacities.capacity / (pv_capacities.area_m2 / 1e6)
    # pv_capacities.loc[pv_capacities.capacity_per_km2 > pv_capacities.capacity_per_km2.quantile(0.95), "capacity_per_km2"] = pv_capacities.capacity_per_km2.quantile(0.95)
    
    fig, ax = plt.subplots()
    gpd.GeoDataFrame(pv_capacities, geometry="geometry").plot(
        column="capacity_per_km2", legend=True, cmap="summer", ax=ax,
        legend_kwds={"label": "Installed Capacity per $km^2$ in MW"}) 
    fig.savefig(filepath.joinpath("pv_de_per_km2.png"))
    
    fig, ax = plt.subplots()
    gpd.GeoDataFrame(pv_capacities, geometry="geometry").plot(
        column="capacity", legend=True, cmap="summer", ax=ax,
        legend_kwds={"label": "Installed Capacity in MW"}) 
    fig.savefig(filepath.joinpath("pv_de.png"))


    
    

    installed_capacities = anymod_installed_capacities(wdir, 2030)
    