
import os
from pathlib import Path
import pomato_data
import pandas as pd 

if __name__ == "__main__":  
    
    os.chdir(Path(pomato_data.__path__[0]).parent )
    settings = {
        "grid_zones": ["DE", "FR", "BE", "LU", "NL"],
        "weather_year": 2019,
        "capacity_year": 2020, 
        # "capacity_year": 2020, 
        "co2_price": 50,
        "split_lines": True,
        # "time_horizon": "01.11.2019 - 30.11.2019",
        "time_horizon": "01.01.2019 - 31.12.2019",
        }
    
    if Path(os.path.abspath("")).name != "pomato_data":
        raise FileNotFoundError("Please Execute the script in the repository itself, use os.chdir() to change path")
    else: 
        wdir = Path(os.path.abspath(""))
    
    data = pomato_data.PomatoData(wdir, settings)
    
    data.inflows.to_csv(wdir.joinpath("inflows.csv"))


    # %%

    timesteps = data.demand_el[["utc_timestamp"]].copy()
    timesteps[["year","week","day"]] = timesteps.utc_timestamp.dt.isocalendar() 
    
    nth_week = 7
    weeks = [2 + i*nth_week for i in range(0, int(timesteps.week.max()/nth_week) + 1)]
    
    timesteps = timesteps[(timesteps.week.isin(weeks))&(timesteps.year == settings["weather_year"])]
    
    # data.storage_level = data.storage_level[data.storage_level.timestep.isin(timesteps.index)]
    
    data.availability = data.availability[data.availability.index.isin(timesteps.index)]
    data.demand_el = data.demand_el[data.demand_el.index.isin(timesteps.index)]
    data.inflows = data.inflows[data.inflows.index.isin(timesteps.index)]

    # t = data.storage_level[data.storage_level.plant == "H1"]
    
    # %% CWE Processing
    data.add_dcline("nNO", "nSE", 3800)
    data.add_dcline("nDK", "nSE", 2000)
    data.add_dcline("nDK", "nNO", 2000)
    data.add_dcline("nCH", "nIT", 2000)
    data.create_basic_ntcs()
    
    data.dclines
    
    # Decommissioning (manual)
    if settings["capacity_year"] == 2030:
        condition_lignite = data.plants.fuel == "lignite"
        condition_coal = data.plants.fuel == "hard coal"
        condition_nuclear = data.plants.fuel == "uran"
        condition_nuclear = data.plants.fuel == "uran"
        condition_gas= data.plants.fuel == "gas"
        condition_de = data.plants.zone == "DE"
        
        data.plants = data.plants.loc[~(condition_lignite & condition_de)]
        data.plants = data.plants.loc[~(condition_nuclear & condition_de)]
        data.plants = data.plants.loc[~(condition_nuclear & condition_de)]
        data.plants = data.plants.loc[~(condition_nuclear & condition_de)]
        data.plants.loc[(condition_gas & condition_de), "g_max"] *= 1.5
        
        data.plants.loc[(condition_nuclear), "g_max"] *= 0.7
        data.plants.loc[(condition_coal), "g_max"] *= 0.7
        data.plants.loc[(condition_lignite), "g_max"] *= 0.7
    
    
    # Remove small plants below a certain threshold 
    if settings["capacity_year"] == 2030:
        threshold = 11
    else:
        threshold = 7
    data.plants[data.plants.g_max > threshold].g_max.sum() / data.plants.g_max.sum()  
    len(data.plants[data.plants.g_max > threshold]) / len(data.plants[data.plants.g_max > 0])
    
    data.plants = data.plants[data.plants.g_max > threshold]
    drop_plants = [p for p in data.availability.columns if p not in data.plants.index]
    data.availability = data.availability.drop(drop_plants, axis=1)
    
    foldername = f"CWE_{settings['capacity_year']}_8_weeks"
    data.save_to_csv(foldername)

    
    # %% Testing 

    data.lines

    # dclines = data.dclines
    # lines = data.lines
    # nodes = data.nodes
    # plants = data.plants
    # zones = data.zones
    # ntc = data.ntc
    # technology = data.technology
    
    # plants_2020 = data.plants.copy()
    # # plants_2030 = data.plants.copy()
    # t = plants_2020[["fuel", "technology", "zone", "g_max"]].groupby(["fuel", "technology", "zone"]).sum().reset_index().fillna(0)
    # t = t.pivot(index="zone", columns=("fuel", "technology"), values="g_max")
    # t.plot.bar(stacked=True)

