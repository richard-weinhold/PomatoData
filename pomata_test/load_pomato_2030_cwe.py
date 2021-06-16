"""CWE Test Case."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pomato 

# Init POMATO with the options file and the dataset
wdir = Path("C:/Users/riw/tubCloud/Uni/Market_Tool/pomato_studies")
mato = pomato.POMATO(wdir=wdir, options_file="profiles/de_2030.json")

mato.load_data(r"C:\Users\riw\Documents\repositories\pomato_data\pomato_datasets\CWE_2030.zip")

mato.options["timeseries"]["redispatch_horizon"] = 168
mato.options["model_horizon"] = [168, 2*168]
mato.options["redispatch"]["include"] = False
mato.options["redispatch"]["zones"] = ["DE", "LU", "FR", "BE", "NL"]
mato.options["fbmc"]["flowbased_region"] = ["DE", "LU", "FR", "BE", "NL"]
mato.options["redispatch"]["cost"] = 50
mato.options["infeasibility"]["electricity"]["cost"] = 500
mato.options["curtailment"]["cost"] = 100
mato.options["grid"]["sensitivity"] = 1e-1
mato.options["plant_types"]["ts"] = ["wind onshore", "wind offshore", "solar"]
mato.options["plant_types"]["es"] = ["hydro_psp"]

mato.data.results = {}
mato.data.missing_data
mato.data.data_validation_report

# Acess the data pre-marketmodel.
nodes = mato.data.nodes
lines = mato.grid.lines
dclines = mato.data.dclines
demand = mato.data.demand_el
zones = mato.data.zones
plants = mato.data.plants

t = plants[plants.plant_type.isna()]

availability = mato.data.availability
ntc = mato.data.ntc

mato._join_julia_instances()

# %% Run the Market Model, including (costbased) Re-Dispatch.
# The Market Result is determined based on the option file.
# The Redispatrch is done to N-0 per default.
# mato.options["type"] = "nodal"
# mato.options["title"] = "Basecase"
# mato.create_grid_representation()
# mato.update_market_model_data()
# mato.run_market_model()

# basecase = mato.data.return_results()
# mato.options["title"] = "FBMC"
# mato.options["fbmc"]["minram"] = 0.4
# mato.options["fbmc"]["lodf_sensitivity"] = 0.15
# mato.options["fbmc"]["cne_sensitivity"] = 0.2
# fb_parameters = mato.create_flowbased_parameters(basecase, enforce_ntc_domain=False)

# mato.options["redispatch"]["include"] = True
# mato.create_grid_representation(flowbased_paramters=fb_parameters)
# mato.update_market_model_data()
# mato.run_market_model()

# mato.options["type"] = "ntc"
# mato.options["title"] = "NTC"
# mato.create_grid_representation()
# mato.update_market_model_data()
# mato.run_market_model()

# %%

result_dir = wdir.joinpath("data_temp/julia_files/results")
folders = [
    result_dir.joinpath("2505_2347_55539"),
    result_dir.joinpath("2505_2356_16086_market_results"),
    result_dir.joinpath("2505_2356_16086_redispatch_DE_LU_FR_BE_NL"),
    result_dir.joinpath("2605_0004_31573_market_results"),
    result_dir.joinpath("2605_0004_31573_redispatch_DE_LU_FR_BE_NL"),
    ]
mato.initialize_market_results(folders)

# %%

market_result = mato.data.results["2505_2356_16086_market_results"]
redispatch_result = mato.data.results["2505_2356_16086_redispatch_DE_LU_FR_BE_NL"]

# gen = redispatch_result.redispatch()
# gen.delta_abs.sum()

# mato.visualization.create_installed_capacity_plot()

# mato.visualization.create_available_intermittent_capacity_plot(mato.data, 
#                                                                zones=["NO"])
# mato.visualization.create_generation_pie(redispatch_result)
# mato.visualization.create_storage_plot(market_result)


#%%
import pomato 
from pomato.visualization import Dashboard
pomato_dashboard = Dashboard(mato)
options = {"debug": False, "use_reloader": False, "port": "8050", "host": "0.0.0.0"}
pomato_dashboard.app.run_server(**options)