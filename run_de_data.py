
#%% Import packages
import os
import sys
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pomato_data import PomatoData

#%%
if __name__ == "__main__":  

    settings = {
        # "grid_zones": ["DE", "FR", "BE", "LU", "NL"],
        "grid_zones": ["DE"],
        "weather_year": 2019,
        # "capacity_year": 2030, 
        "capacity_year": 2020, 
        "co2_price": 60,
        "split_lines": True,
        # "time_horizon": "01.11.2019 - 30.11.2019",
        "time_horizon": "01.05.2019 - 31.05.2019",
        }
    
    wdir = Path(os.path.dirname(os.path.abspath(__file__)))
    data = PomatoData(wdir, settings)

    # %% DE Processing
    data.add_dcline("nDK", "nSE", 2000)
    data.create_basic_ntcs()

    tr = 13

    data.plants[data.plants.g_max > tr].g_max.sum() / data.plants.g_max.sum()  
    len(data.plants[data.plants.g_max > tr]) / len(data.plants[data.plants.g_max > 0])
    
    data.plants = data.plants[data.plants.g_max > tr]
    drop_plants = [p for p in data.availability.columns if p not in data.plants.index]
    data.availability = data.availability.drop(drop_plants, axis=1)
    
    foldername = f"DE_{settings['capacity_year']}"
    data.save_to_csv(foldername)
    
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
    
    # plants_2020 = data.plants.copy()
    # # plants_2030 = data.plants.copy()
    # t = plants_2020[["fuel", "technology", "zone", "g_max"]].groupby(["fuel", "technology", "zone"]).sum().reset_index().fillna(0)
    # t = t.pivot(index="zone", columns=("fuel", "technology"), values="g_max")
    # t.plot.bar(stacked=True)


# %%
