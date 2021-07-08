
import os
from pathlib import Path
import pomato_data

if __name__ == "__main__":  

    settings = {
        "grid_zones": ["DE", "FR", "BE", "LU", "NL"],
        "weather_year": 2019,
        # "capacity_year": 2030, 
        "capacity_year": 2020, 
        "co2_price": 60,
        "split_lines": False,
        # "time_horizon": "01.11.2019 - 30.11.2019",
        "time_horizon": "01.03.2019 - 02.03.2019",
        }
    
    if Path(os.path.abspath("")).name != "pomato_data":
        raise FileNotFoundError("Please Execute the script in the repository itself, use os.chdir() to change path")
    else: 
        wdir = Path(os.path.abspath(""))
    
    data = pomato_data.PomatoData(wdir, settings)

    
    # %% CWE Processing
    data.add_dcline("nNO", "nSE", 3800)
    data.add_dcline("nDK", "nSE", 2000)
    data.add_dcline("nDK", "nNO", 2000)
    data.add_dcline("nCH", "nIT", 2000)
    data.create_basic_ntcs()
    
    data.dclines
    
    # Decommissioning (manual)
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
    threshold = 12
    data.plants[data.plants.g_max > threshold].g_max.sum() / data.plants.g_max.sum()  
    len(data.plants[data.plants.g_max > threshold]) / len(data.plants[data.plants.g_max > 0])
    
    data.plants = data.plants[data.plants.g_max > threshold]
    drop_plants = [p for p in data.availability.columns if p not in data.plants.index]
    data.availability = data.availability.drop(drop_plants, axis=1)
    
    foldername = f"CWE_{settings['capacity_year']}_nov_decarb"
    data.save_to_csv(foldername)
    
    # %% Testing 
    # availability = data.availability
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

