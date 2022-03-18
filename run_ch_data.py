# %%
import os
import pandas as pd
from pathlib import Path
from pomato_data.pomato_data import PomatoData
import pomato_data

if __name__ == "__main__":  
    
    os.chdir(Path(pomato_data.__path__[0]).parent)
    if Path(os.path.abspath("")).name != "pomato_data":
        raise FileNotFoundError("Please Execute the script in the repository itself, use os.chdir() to change path")
    else: 
        wdir = Path(os.path.abspath(""))
        
    settings = {
        "grid_zones": ["CH", "AT", "DE", "IT", "FR"],
        "except_zones": ["MT", "ME"],
        "include_neighbors": False,
        # "grid_zones": ["CH", "AT"],
        "weather_year": 2019,
        "capacity_year": 2035, 
        "co2_price": 100,
        "split_lines": False,
        "capacity_source": "manual",
        "capacity_file": "data_in/res/manual_capacities_CH_2035.csv",
        "time_horizon": "01.01.2019 - 31.12.2019",
        }
    
    data = PomatoData(wdir, settings)
    
# %%
    
    # # %% DE Processing
    # data.add_dcline("nDK", "nSE", 2000)
    # data.create_basic_ntcs()
    # data.inflows["H1"].plot()
    
    # # add heat
    data.plants[["h_max", "zone"]].groupby("zone").sum()
    # heat = pd.read_csv("data_in/heat/reference_data_set-v1.0.0.csv", skiprows=6, index_col=0, thousands=",").astype(float)
    # heat["DE"] = heat.sum(axis=1)
    data.demand_h = pd.DataFrame()
    # data.demand_h["timestep"] = data.demand_el.index
    # data.demand_h["heatarea"] = "DE"
    # data.demand_h["demand_h"] = heat.loc[data.demand_el.index, "DE"].values
    data.heatareas = pd.DataFrame()
    # data.plants.loc[(data.plants.zone == "DE")&(data.plants.h_max >0), "heatarea"] = "DE"
    

    # %% Demand Scaling EL
    
    filepath = wdir.joinpath("data_in/demand/manual_demand_2035.csv")
    zonal_demand_2035 = pd.read_csv(filepath, index_col=0)
    zonal_demand_2035.loc["CH"] = 67.4e3 # value from energieperspectiven 50+
    
    for zone in data.zones.index:
        tmp_nodes = data.nodes[data.nodes.zone == zone].index
        zonal_demand = data.demand_el.loc[:, tmp_nodes]
        demand_scaling = zonal_demand_2035.loc[zone, "value"]*1e3/(zonal_demand.sum().sum())
        data.demand_el.loc[:, tmp_nodes] *= demand_scaling
        
    # check
    de_nodes = data.nodes[data.nodes.zone == "DE"].index
    de_demand = data.demand_el.loc[:, de_nodes]
    de_demand.sum().sum()/1e6
    
    # # Demand scaling heat
    cond_de = data.plants.zone == "DE"
    cond_ch = data.plants.zone == "CH"
    cond_nuke = data.plants.fuel == "uran"
    cond_biomass = data.plants.fuel == "biomass"
    cond_biomass = data.plants.fuel == "biomass"
    cond_ror = data.plants.technology == "ror"
    
    data.plants["availability"] = 1
    data.plants.loc[cond_de & cond_biomass, "availability"] = 0.68
    data.plants.loc[cond_ch & cond_biomass, "availability"] = 0.68
    data.plants.loc[cond_de & cond_nuke, "availability"] = 0.95
    data.plants.loc[cond_ch & cond_nuke, "availability"] = 0.95
    data.plants.loc[cond_ch & cond_ror, "availability"] = 0.52
    
    # %% Remove small plants below a certain threshold 
    threshold = 10
    data.plants[data.plants.g_max > threshold].g_max.sum() / data.plants.g_max.sum()  
    len(data.plants[data.plants.g_max > threshold]) / len(data.plants[data.plants.g_max > 0])
    
    data.plants = data.plants[data.plants.g_max > threshold]
    drop_plants = [p for p in data.availability.columns if p not in data.plants.index]
    data.availability = data.availability.drop(drop_plants, axis=1)

    # %% Hydro 
    data.plants.loc[:, "storage_capacity"] *= 0.8
    data.plants.loc["H23", "storage_capacity"] = 225220
    data.plants.loc["H1", "storage_capacity"] = 400170
    
    data.plants[data.plants.zone == "CH"].storage_capacity.sum()
    data.plants.plant_type.unique()
    es = ["hydro_psp", "hydro_res"]
    data.plants.loc["H1054"]

    # Psp without d_max > d_max == g_max
    cond_psp = (data.plants.plant_type == "hydro_psp")
    data.plants.loc[cond_psp&(data.plants.d_max==0), "d_max"] = data.plants.loc[cond_psp&data.plants.d_max==0, "g_max"] 
    data.plants.loc[(data.plants.eta > 0.8)&(data.plants.d_max > 0), "eta"] = 0.75

    # Psp less than 12h gen > no inflow
    cond_small_reservoir = data.plants.g_max*168 > data.plants.storage_capacity.fillna(0)
    for p in data.plants[cond_psp&cond_small_reservoir].index:
        if p in data.inflows.columns:
            # data.inflows.loc[:, p] = 0
            data.inflows.drop(p, axis=1, inplace=True)
    
    # Reduce inflow of PSP plants, as its inflated by charging
    for p,x in (data.plants.loc[cond_psp, "d_max"] / 
                ( data.plants.loc[cond_psp, "d_max"] + data.plants.loc[cond_psp, "g_max"])).items():
            if p in data.inflows.columns:
                data.inflows.loc[:, p] *= x 
        
    # %% Highly connected nodes
    
    lines_i = data.lines[data.lines.technology != "transformer"].copy()
    lines_j = data.lines[data.lines.technology != "transformer"].copy()
    lines_i = lines_i[["node_i", "capacity"]].reset_index().rename(columns={"node_i": "node"})
    lines_j = lines_j[["node_j", "capacity"]].reset_index().rename(columns={"node_j": "node"})
    lines = pd.concat([lines_i, lines_j])
    lines_sum = lines.groupby("node").sum()
    lines_sum["zone"] = data.nodes.loc[lines_sum.index, "zone"] 
    lines_sum = lines_sum.sort_values(by="capacity", ascending=False)
    
    # add batteries
    filepath = wdir.joinpath("data_in/res/manual_capacities_CH_2035.csv")

    capacities = pd.read_csv(filepath)
    capacities = capacities[capacities.country.isin(data.zones.index)&(capacities.group=="battery")]
    cols = list(data.plants.columns)
    add_batteries = []
    add_batteries = []
    nodes = data.nodes.copy()
    nodes.loc[nodes.name.isna(), "name"] = ""
    for zone in data.zones.index:
        tmp_nodes = lines_sum[lines_sum.zone == zone]
        battery_cap = sum(capacities.loc[capacities.country == zone, "value"])
        if battery_cap == 0:
            break
        
        number = 10 # top ten connected nodes
        for n in tmp_nodes.index[:number]:
            tmp_plant = {
                "mc_el": 10,
                "mc_heat": 0,
                "eta": 0.8,
                "node": n,
                "g_max": battery_cap/number,
                "d_max": battery_cap/number,
                "h_max": 0,
                "status": "online",
                "availability": 1,
                "zone": data.nodes.loc[n, "zone"],
                "plant_type": "battery",
                "technology": "battery",
                "fuel": "electricity",
                "storage_capacity": battery_cap/number*24,
                "lat": data.nodes.loc[n, "lat"],
                "lon": data.nodes.loc[n, "lon"],
                "name": "battery_" + nodes.loc[n, "name"],
                "commissioned": 2020
            }
            add_batteries.append([tmp_plant[d] if d in tmp_plant.keys() else None for d in cols])
    
    batteries = pd.DataFrame(columns=cols, data=add_batteries)
    batteries = batteries.set_index("battery_" + batteries.node)
    data.plants = pd.concat([data.plants, batteries])
    
    # %% decomm
    
    # decomm = pd.read_csv(wdir.joinpath("data_in/plants/decomm_2022.csv"), index_col=0, encoding="ISO-8859-1")
    # data.plants = data.plants.loc[~data.plants.index.isin(decomm.index)]    
    
    # %%
    from pomato_data.res.capacities import other_res_capacities
    from pomato_data.auxiliary import get_countries_regions_ffe
    import geopandas as gpd
    import shapely
    # add 4GW other res for non-kkw
    # add 400 MW other res for kkw
    other_res_capacity = 4000
    number_of_plants = 5

    # if other_res_capacity == 400:
        # CH KKWs
    data.plants[(data.plants.zone == "CH")&(data.plants.fuel == "uran")]
    data.plants = data.plants.drop(['p1719', 'p1721', 'p1722'])
        
    # use current installations as distribution
    other_res = other_res_capacities(wdir)
    ch_other_res = other_res[other_res.nuts_id.str.contains("CH")&(other_res.fuel != "hydro")].groupby("nuts_id").sum()
    ch_other_res *= other_res_capacity/ch_other_res.sum()

    country_data, nuts_data = get_countries_regions_ffe()    
    nuts_data = gpd.GeoDataFrame(nuts_data[nuts_data.country == "CH"]).set_crs("EPSG:4326")
    
    nodes = data.nodes.copy()
    geometry = [shapely.geometry.Point(xy) for xy in zip(nodes.lon, nodes.lat)]
    condition = [n.within(country_data.loc[zone, "geometry"]) for n in geometry]
    nodes = gpd.GeoDataFrame(nodes, crs="EPSG:4326", geometry=geometry).loc[condition]
    
    capacity_nodes = nodes[['voltage', 'name', 'zone', "info", "lat", "lon"]].copy()
    nuts_to_nodes = gpd.sjoin(nodes, nuts_data, how='left', predicate='within')
    
    new_other_res_plants = pd.DataFrame()
    for nuts in ch_other_res.index:
        tmp_nodes = lines_sum[lines_sum.index.isin(nuts_to_nodes[nuts_to_nodes.name_short == nuts].index)].index[:number_of_plants]
        tmp_other_res_plants = capacity_nodes.loc[tmp_nodes, ["zone", 'name', "lat", "lon"]].reset_index()
        tmp_other_res_plants["index"] = tmp_other_res_plants.node.astype(str) + "_" + "other_res"
        tmp_other_res_plants["name"] = tmp_other_res_plants.node.astype(str) + "_" + "other_res"
        tmp_other_res_plants = tmp_other_res_plants.set_index("index")
        tmp_other_res_plants[["fuel", "technology"]] = "biomass", "ccgt"
        tmp_other_res_plants["eta"] = data.technology.set_index(["fuel", "technology"]).loc[("biomass", "ccgt"), "eta"]
        tmp_other_res_plants["plant_type"] = "other_res"
        tmp_other_res_plants["commissioned"] = 1900
        tmp_other_res_plants["status"] = "online"
        tmp_other_res_plants[["chp", "heatarea", "city", "postcode", "street"]] = None
        tmp_other_res_plants["g_max"] = ch_other_res.loc[nuts, "capacity"] / number_of_plants
        new_other_res_plants = pd.concat([new_other_res_plants, tmp_other_res_plants])  
    new_other_res_plants.loc[:, "g_max"] *= other_res_capacity / new_other_res_plants.loc[:, "g_max"].sum()
    data.plants = pd.concat([data.plants, new_other_res_plants])

    # %%
    
    # %% Two Gas plants 
    # CH,gas,500.0
    number_of_plants = 3
    
    gas_nodes = lines_sum[lines_sum.zone == "CH"].index[:number_of_plants]
    tmp_gas_plants = capacity_nodes.loc[gas_nodes, ["zone", 'name', "lat", "lon"]].reset_index()
    tmp_gas_plants["index"] = tmp_gas_plants.node.astype(str) + "_" + "gas"
    tmp_gas_plants["name"] = tmp_gas_plants.node.astype(str) + "_" + "gas"
    tmp_gas_plants = tmp_gas_plants.set_index("index")
    tmp_gas_plants[["fuel", "technology"]] = "gas", "ccgt"
    tmp_gas_plants["eta"] = data.technology.set_index(["fuel", "technology"]).loc[("gas", "ccgt"), "eta"]
    tmp_gas_plants["plant_type"] = "conventional"
    tmp_gas_plants["commissioned"] = 1900
    tmp_gas_plants["status"] = "online"
    tmp_gas_plants[["chp", "heatarea", "city", "postcode", "street"]] = None
    tmp_gas_plants["g_max"] = 500 / number_of_plants
    data.plants = pd.concat([data.plants, tmp_gas_plants])
    data.marginal_costs()
    data.uniquify_marginal_costs()

    data.plants[data.plants.mc_el.isna()]
    
    # %%
    
    plants = data.plants.copy()
    t = plants[["fuel", "technology", "zone", "g_max"]].groupby(["fuel", "technology", "zone"]).sum().reset_index().fillna(0)
    t = t.pivot(index="zone", columns=("fuel", "technology"), values="g_max")
    t.plot.bar(stacked=True)
    t.loc["CH"]
    
    # %% Exchnage: First values from TYNDP/AnyMOD
    
    filename = "data_in/anymod_results/results_exchange_base_120_full_202201271632.csv"
    anymod_exchnage = pd.read_csv(wdir.joinpath(filename))
    
    anymod_exchnage = anymod_exchnage[anymod_exchnage.variable == "capaExc"]
    anymod_exchnage = anymod_exchnage[anymod_exchnage.exchange.str.contains("powerGrid")]
    anymod_exchnage.loc[:, "zone_from"] = anymod_exchnage["region_from"].str.split("<", expand=True)[1].values
    anymod_exchnage.loc[:, "zone_to"] = anymod_exchnage["region_to"].str.split("<", expand=True)[1].values
    anymod_exchnage.loc[:, "zone_from"] = anymod_exchnage["zone_from"].str.strip(" ")
    anymod_exchnage.loc[:, "zone_to"] = anymod_exchnage["zone_to"].str.strip(" ")
    
    year = '2030'
    anymod_exchnage = anymod_exchnage[anymod_exchnage.timestep_superordinate_dispatch.str.contains(str(year))]
    anymod_exchnage.loc[:, "value"] = anymod_exchnage["value"]*1000
    
    anymod_exchnage = anymod_exchnage[~(anymod_exchnage.zone_from == anymod_exchnage.zone_to)]
    anymod_exchnage = anymod_exchnage[["zone_from", "zone_to", "value"]].groupby(["zone_from", "zone_to"]).sum()

    for row in data.ntc.iterrows():
        if (row[1].zone_i, row[1].zone_j) in anymod_exchnage.index:
            data.ntc.loc[row[0], "ntc"] = anymod_exchnage.loc[(row[1].zone_i, row[1].zone_j), "value"]
    
    # Second Values from FE Stromzusammenarbeit
    filename = "data_in/exchange/ch_ntcs.csv"

    ch_ntc = pd.read_csv(wdir.joinpath(filename))
    ch_ntc = ch_ntc[ch_ntc.scenario == "no-cooperation"].set_index(["zone_i", "zone_j"])
    for row in data.ntc.iterrows():
        if (row[1].zone_i, row[1].zone_j) in ch_ntc.index:
            data.ntc.loc[row[0], "ntc"] = ch_ntc.loc[(row[1].zone_i, row[1].zone_j), "ntc"]
            
    data.set_ntc_for_no_connection(set_min_values=False)
    
    
    # %%
    condition_lignite = data.plants.fuel == "lignite"
    condition_coal = data.plants.fuel == "hard coal"
    condition_nuclear = data.plants.fuel == "uran"
    condition_gas= data.plants.fuel == "gas"
    condition_de = data.plants.zone == "DE"
    data.plants = data.plants.loc[~(condition_lignite & condition_de)]
    data.plants = data.plants.loc[~(condition_nuclear & condition_de)]
    data.plants = data.plants.loc[~(condition_coal & condition_de)]
    
    # Sum GAS DE 33GW 
    sum_gas = data.plants.loc[(condition_gas & condition_de), "g_max"].sum()    
    data.plants.loc[(condition_gas & condition_de), "g_max"] *= 33e3/sum_gas
    
    # Sum GAS IT 47GW
    condition_it = data.plants.zone == "IT"
    sum_gas = data.plants.loc[(condition_gas & condition_it), "g_max"].sum()    
    data.plants.loc[(condition_gas & condition_it), "g_max"] *= 47e3/sum_gas
    
    # Sum Nuclear FR 47GW
    condition_fr = data.plants.zone == "FR"
    sum_kkw = data.plants.loc[(condition_nuclear & condition_fr), "g_max"].sum()    
    data.plants.loc[(condition_nuclear & condition_fr), "g_max"] *= 57e3/sum_kkw

    
    # %%
    
    foldername = f"CH_{settings['capacity_year']}_v04"
    data.save_to_csv(foldername, 
                     path=Path(r"C:\Users\riw\tubCloud\Uni\Market_Tool\pomato_studies\data_input"))
    
    # %% Testing 
    # # availability = data.availability
    # demand = data.demand_el
    
    # dclines = data.dclines
    # lines = data.lines
    # nodes = data.nodes
    # plants = data.plants
    # zones = data.zones
    # ntc = data.ntc
    # technology = data.technology
    # tt = data.plants.loc[cond_de & cond_nuke, :]
        
    # data.plants.columns
    
    # plants_2020 = data.plants.copy()
    # # plants_2030 = data.plants.copy()
    # t = plants_2020[["fuel", "technology", "zone", "g_max"]].groupby(["fuel", "technology", "zone"]).sum().reset_index().fillna(0)
    # t = plants_2020[["fuel", "zone", "g_max"]].groupby(["fuel", "zone"]).sum().reset_index().fillna(0)
    # t = t.pivot(index="zone", columns=("fuel"), values="g_max")
    # t = t.pivot(index="zone", columns=("fuel", "technology"), values="g_max")
    # t.plot.bar(stacked=True)
    # t.loc["DE"]
