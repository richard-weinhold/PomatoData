import pandas as pd
import numpy as np
from pathlib import Path
### Renaming Dictionaries
# convert opsd fuels to elmod fuels
fuel_dict = {"natural gas": "gas", "biomass and biogas": "biomass", 'bioenergy': "biomass",
             "other bioenergy and renewable waste": "biomass",
             "nuclear": "uran",
             "other fuels": "waste", "other fossil fuels": "waste", "fossil fuels": "gas",
             "other fossil fuels": "waste", "non-renewable waste": "waste",
             "mixed fossil fuels": "waste", "other or unspecified energy sources": "waste",
             'solar': "sun", 'geothermal': "sun"}

# convert opsd technologies to elmod technologies
tech_dict = {"steam turbine": "steam", "waste": "steam", "other fossil fuels": "steam",
             "combustion engine": "generator",
             "combined cycle": "ccgt", "ccgt and biogas": "ccgt", "biomass": "ccgt",
             "sewage and landfill gas": "ccgt", "biomass and biogas": "ccgt",
             "gas turbine": "ccgt",
             "run-of-river": "ror",
             "pumped storage": "psp", "storage": "psp", "storage technologies": "psp",
             "pumped storage with natural inflow": "psp",
             "reservoir technologies": "reservoir",
             "hydro": "ror",
             "photovoltaics": "solar", "photovoltaics ground": "solar",
             "offshore": "wind offshore", "wind_offshore": "wind offshore",
             "wind": "wind onshore", "onshore": "wind onshore", "wind_onshore": "wind onshore",
             }

def plants_naming(plants):
    """Adjust names so that tech/fuel data from dynElmod and plant data from OPSD are coherent"""
    # convert to lower case for processable contents
    for col in ["fuel", "technology"]:
        plants[col] = plants[col].str.lower()
    plants.fuel = plants.fuel.str.lstrip(" ").str.rstrip(" ")
    plants.technology = plants.technology.str.lstrip(" ").str.rstrip(" ")

    ## Renaming Fuels and technologies with Dics at top of code
    plants.technology = plants.technology.replace(tech_dict, regex=False)
    plants.fuel = plants.fuel.replace(fuel_dict, regex=False)
    
    # fix things
    plants.loc[plants.name.str.contains("COO").fillna(False), "technology"] = "psp"

    plants.loc[(plants.technology.isnull())&(plants.fuel == "hydro"), "technology"] = "ror"
    plants.loc[(plants.technology.isnull())&(plants.fuel == "sun"), "technology"] = "solar"
    plants.loc[(plants.technology.isnull())&(plants.fuel == "wind"), "technology"] = "wind onshore"
    plants.loc[(plants.technology.isnull())&(plants.fuel == "biomass"), "technology"] = "ccgt"
    plants.loc[(plants.technology.isnull())&(plants.fuel == "gas"), "technology"] = "ccgt"
    plants.loc[(plants.technology.isnull())&(plants.fuel == "oil"), "technology"] = "ccot"
    
    # Rest is steam    
    plants.loc[plants.technology.isnull(), "technology"] = "steam"
    plants.loc[plants.fuel == "waste", "technology"].unique()
    
    plants.loc[plants.fuel == "hard coal", "technology"] = "steam"
    plants.loc[plants.fuel == "waste", "technology"] = "steam"
    plants.loc[plants.fuel == "oil", "technology"] = "ccot"
    plants.loc[plants.technology == "psp", "fuel"] = "hydro"
    plants.loc[plants.fuel ==  "storage", "technology"] = "psp"
    plants.loc[plants.fuel ==  "storage", "fuel"] = "hydro"

    return plants

def efficiency_estimate(filepath, plants):
    """ Replace efficiencies with estimates fram literature"""
    file_path = str(filepath.joinpath("data_in/plants/input_efficiency_literature_by_fuel_technology.csv"))
    efficiency_data = pd.read_csv(file_path, header=0, sep=",")
    for col in ["technology", "fuel"]:
        efficiency_data[col] = efficiency_data[col].str.lower()
    # Rename
    efficiency_data.technology = efficiency_data.technology.replace(tech_dict, regex=False)
    efficiency_data.fuel = efficiency_data.fuel.replace(fuel_dict, regex=False)

    efficiency_data = efficiency_data.groupby(["technology", "fuel"]).mean().reset_index()
    tmp = plants[["technology", "fuel", "commissioned"]]
    tmp = pd.merge(tmp, efficiency_data, on=["technology", "fuel"], how="left").set_index(tmp.index)
    tmp["eta_estimate"] = tmp.efficiency_intercept + tmp.commissioned*tmp.efficiency_slope
    plants.loc[tmp.index[~tmp.eta_estimate.isna()], "eta"] = tmp.eta_estimate[~tmp.eta_estimate.isna()]
    return plants

def efficiency(plants, technology, fuel):
    """Calculate the efficiency for plants that dont have it maunually set.

    Calculates/assigns the efficiency of a power plant based on information in tech/fuel/plant
    tables.

    If none is available eta defaults to the default value set in the option file.
    """
    tmp = plants[['technology', 'fuel', 'plant_type']][plants.eta.isnull()]
    tmp = pd.merge(tmp, technology[['technology', 'fuel', 'eta']],
                   how='left', on=['technology', 'fuel']).set_index(tmp.index)
    plants.loc[plants.eta.isnull(), "eta"] = tmp.eta
    # ## If there are still data.plants without efficiency
    # default_value = 0.3
    # plants.loc[plants.eta.isnull(), "eta"] = default_value
    return plants 

def process_plants(filepath):
    
    plants = pd.read_csv(filepath.joinpath("data_in/plants/plants.csv"), index_col=0)
    plants.zone = plants.zone.str.lstrip(" ").str.rstrip(" ")
    plants = plants.rename(columns={"capacity": "g_max", "chp_capacity": "h_max", "country": "zone"})

    fuel = pd.read_csv(filepath.joinpath("data_out/fuel/fuel.csv"), index_col=0)
    technology = pd.read_csv(filepath.joinpath("data_out/technology/technology.csv"), index_col=0)
    
    plants = plants_naming(plants)
    
    plants.loc[:, "plant_type"] = pd.merge(plants[["fuel", "technology"]], technology, 
                                           how="left", on=["fuel", "technology"])["plant_type"].values    
    
    plants = efficiency_estimate(filepath, plants)
    plants = efficiency(plants, technology, fuel)
    return plants 

# %%
if __name__ == "__main":
    
    filepath = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    plants = process_plants(filepath)
    plants.to_csv(filepath.joinpath('data_out/plants/plants.csv'))
    df = plants[plants.eta.isna()]

    

