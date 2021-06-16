# -*- coding: utf-8 -*-
"""
Created on Wed May 26 15:19:26 2021

@author: s792
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pomato import POMATO
import plotly.graph_objects as go
import plotly
import plotly.express as px
import os
# from plots import *
# from auxiliary import *
from more_itertools import distinct_combinations
from pomato.visualization import Dashboard, Visualization
import copy

import plotly.io as pio

pio.renderers.default='browser'

#%% Create mato object
wdir = Path(r"C:\Users\s792\Documents\pomato_data\pomata_test")
mato = POMATO(wdir=wdir, options_file="data_input/pomato_2030.json")
mato.load_data('data_input/CWE_2030/')

#%% Read in data
mato.initialize_market_results([Path(r"C:\Users\s792\Documents\pomato_data\pomata_test\data_temp\julia_files\results\2605_1459_04455_market_results"),
                                Path(r"C:\Users\s792\Documents\pomato_data\pomata_test\data_temp\julia_files\results\2605_1459_04455_redispatch_LU_ES_SE_AT_BE_CZ_DK_FR_DE_IT_NL_NO_PL_CH_UK")])

#%% Plotting
res_names = list(mato.data.results.keys())
mato.visualization.create_generation_plot(mato.data.results[res_names[0]], nodes=None)

mato.start_dashboard(debug=True)

pomato_dashboard = Dashboard(mato)
app = pomato_dashboard.app
server = pomato_dashboard.app.server
options = {"debug": True, "use_reloader": False, "port": "8050", "host": "0.0.0.0"}
app.run_server(**options)
