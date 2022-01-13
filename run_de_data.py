
import os
import pandas as pd
from pathlib import Path
from pomato_data.pomato_data import PomatoData
import pomato_data

if __name__ == "__main__":  
    
    os.chdir(Path(pomato_data.__path__[0]).parent )
    if Path(os.path.abspath("")).name != "pomato_data":
        raise FileNotFoundError("Please Execute the script in the repository itself, use os.chdir() to change path")
    else: 
        wdir = Path(os.path.abspath(""))
        
    settings = {
        "grid_zones": ["DE"],
        "weather_year": 2019,
        "capacity_year": 2022, 
        "co2_price": 60,
        "split_lines": False,
        "capacity_source": "manual",
        "capacity_file": wdir.joinpath("data_in/res/manual_capacities_2022.csv"),
        "time_horizon": "01.01.2019 - 31.12.2019",
        }
    
        
    data = PomatoData(wdir, settings)
    
    # data.ntc[data.ntc.zone_i == "PL"]
    data.lines.voltage.unique()
    data.lines["count"] = 1
    data.lines.groupby("voltage").sum()["count"]
    
    data.plants.plant_type.unique()
    data.plants["count"] = 1
data.nodes
    data.plants[(data.plants.plant_type == "conventional")&(data.plants.zone=="DE")].g_max.sum()
    data.plants[(data.plants.zone=="DE")].groupby("plant_type").sum()[["count", "g_max"]]
    
    # %% DE Processing
    data.add_dcline("nDK", "nSE", 2000)
    data.create_basic_ntcs()
    
    # add heat
    heat = pd.read_csv("data_in/heat/reference_data_set-v1.0.0.csv", skiprows=6, index_col=0, thousands=",").astype(float)
    heat["DE"] = heat.sum(axis=1)
    data.demand_h = pd.DataFrame()
    data.demand_h["timestep"] = data.demand_el.index
    data.demand_h["heatarea"] = "DE"
    data.demand_h["demand_h"] = heat.loc[data.demand_el.index, "DE"].values
    data.heatareas = pd.DataFrame(index=["DE"])
    data.plants.loc[(data.plants.zone == "DE")&(data.plants.h_max >0), "heatarea"] = "DE"
    
    # %%

    # Demand Scaling EL
    de_nodes = data.nodes[data.nodes.zone == "DE"].index
    demand_de = data.demand_el.loc[:, de_nodes]
    demand_scaling = 570.8/(demand_de.sum().sum()/1e6)
    data.demand_el.loc[:, de_nodes] *= demand_scaling
    
    # Demand scaling heat
    cond_de = data.plants.zone == "DE"
    cond_nuke = data.plants.fuel == "uran"
    cond_biomass = data.plants.fuel == "biomass"
    
    data.plants["availability"] = 1
    data.plants.loc[cond_de & cond_biomass, "availability"] = 0.68
    data.plants.loc[cond_de & cond_nuke, "availability"] = 0.95
    
    # %% Remove small plants below a certain threshold 
    threshold = 10
    data.plants[data.plants.g_max > threshold].g_max.sum() / data.plants.g_max.sum()  
    len(data.plants[data.plants.g_max > threshold]) / len(data.plants[data.plants.g_max > 0])
    
    data.plants = data.plants[data.plants.g_max > threshold]
    drop_plants = [p for p in data.availability.columns if p not in data.plants.index]
    data.availability = data.availability.drop(drop_plants, axis=1)
    
    # %% decomm
    
    decomm = pd.read_csv(wdir.joinpath("data_in/plants/decomm_2022.csv"), index_col=0)
    data.plants = data.plants.loc[~data.plants.index.isin(decomm.index)]    
    
    # %%
    # if settings["capacity_year"] == 2030:
    #     # Decommissioning (manual)

    #     condition_lignite = data.plants.fuel == "lignite"
    #     condition_coal = data.plants.fuel == "hard coal"
    #     condition_nuclear = data.plants.fuel == "uran"
    #     condition_gas= data.plants.fuel == "gas"
    #     condition_de = data.plants.zone == "DE"
    #     data.plants = data.plants.loc[~(condition_lignite & condition_de)]
    #     data.plants = data.plants.loc[~(condition_nuclear & condition_de)]
    #     data.plants.loc[(condition_gas & condition_de), "g_max"] *= 1.5
    #     data.plants.loc[(condition_nuclear|condition_coal|condition_lignite), "g_max"] *= 0.7
    
    foldername = f"DE_{settings['capacity_year']}_Wochenbericht_heat"
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
    tt = data.plants.loc[cond_de & cond_nuke, :]
        
    data.plants.columns
    
    plants_2020 = data.plants.copy()
    # # plants_2030 = data.plants.copy()
    t = plants_2020[["fuel", "technology", "zone", "g_max"]].groupby(["fuel", "technology", "zone"]).sum().reset_index().fillna(0)
    # t = plants_2020[["fuel", "zone", "g_max"]].groupby(["fuel", "zone"]).sum().reset_index().fillna(0)
    # t = t.pivot(index="zone", columns=("fuel"), values="g_max")
    t = t.pivot(index="zone", columns=("fuel", "technology"), values="g_max")
    t.plot.bar(stacked=True)
    t.loc["DE"]
