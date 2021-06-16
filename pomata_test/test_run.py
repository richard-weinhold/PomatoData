# -*- coding: utf-8 -*-
"""
Created on Wed May 26 10:25:55 2021

@author: s792
"""

import pandas as pd
from pomato import POMATO
import math
import geopandas as gpd
import os
from pathlib import Path
from auxiliary import *
import plotly.graph_objects as go
import pomato

#%% Read in data
wdir = Path(r"C:\Users\s792\Documents\pomato_data\pomata_test")
mato = POMATO(wdir=wdir, options_file="data_input/pomato_2030.json")
mato.load_data('data_input/CWE_2030/')

#%% Set pomato options
mato.options["type"] = "ntc"
mato.options["model_horizon"] = [704, 720]
mato.options["timeseries"]["redispatch_horizon"] = 16
mato.options["timeseries"]["market_horizon"] = 16
mato.options["redispatch"]["include"] = True
mato.options["redispatch"]["zonal_redispatch"] = False
mato.options["redispatch"]["zones"] = mato.data.zones.index.tolist()
mato.options["infeasibility"]["electricity"]["bound"] = 2000
mato.options["infeasibility"]["electricity"]["cost"] = 1000
mato.options["redispatch"]["cost"] = 50
mato.options["curtailment"]["include"] = True
mato.options["curtailment"]["cost"] = 100
mato.options["grid"]["capacity_multiplier"] = 1

#%% Run 
mato.create_grid_representation()
mato.update_market_model_data()

mato.data.results = dict()
mato.run_market_model()
mato._join_julia_instances()
