<img  height="24" src="https://raw.githubusercontent.com/richard-weinhold/pomato/main/docs/_static/graphics/pomato_logo_small.png"> PomatoData <img  height="24" src="https://raw.githubusercontent.com/richard-weinhold/pomato/main/docs/_static/graphics/pomato_logo_small.png">
=========================================================================================================================================================

Overview
--------

PomatoData provides means to generate a compatible data-set to be used with [POMATO](https://github.com/richard-weinhold/pomato/). PomatoData uses various ressources of the open data community, which are listed below. 

To compile a data-set this repository holds multiple scripts that can be categorized into 1) pre-processing the raw data and 2) processing the final data-set. This is reflected in the folder structure where *data_in* holds raw, unaltered data and *data_out* holds processed data. 

The script *preprocess_input_data.py* includes all preprocessing, including the creation of a cutout of weather data that is used to compute availability and inflow timeseries. This can take substantial amounts of time (hours), but has only to be done once. 

Note, the raw input data has to be aquired by the user manually, sources are included in *preprocess_input_data.py* and the specified location within *data_in*. 

Data Sources
------------

The dataset is compiled from different fantastic projects in the open data community. Namely these are:

   - Conventional power plant data is taken from the [Open Power System Data Platform (OPSD)](https://open-power-system-data.org/).
   - Geographic information is used from the [ExtremOS](https://opendata.ffe.de/project/extremos/) project of Forschungsstelle für Energiewirtschaft (FfE) and their FfE Open Data Portal.   
   - Wind and PV capacities are distributed using the NUTS3 Potentials from FfE. 
   - Future capacities are taken from results of the large scale energy system model [AnyMod](https://github.com/leonardgoeke/AnyMOD.jl)
   - NUTS3 availability timeseries for wind and solar are generated using the [atlite](https://github.com/PyPSA/atlite), package. Offshore availability based on EEZ regions of FfE.
   - The grid data comes from the [GridKit](https://github.com/bdw/GridKit) project, more specifically [PyPSA/pypsa-eur](https://github.com/PyPSA/pypsa-eur/tree/master/data/entsoegridkit) fork, which contains more recent data. 
   - Hydro Plants are taken from the [JRC Hydro-power plants database](https://github.com/energy-modelling-toolkit/hydro-power-database) and inflows are determined using the [atlite](https://github.com/PyPSA/atlite) hydro capabilities and scaled using annual generation.
   - Load, commercial exchange and weekly storage level from [ENTSO-E Transparency platform](https://transparency.entsoe.eu/)


Example
-------
<img  src="https://github.com/richard-weinhold/pomato_data/blob/main/docs/_static/generation_plot.png?raw=true">

After all data is preprocessed, PomatoData can be used to generate data sets. Two scripts *run_cwe_data.py* and *run_de_data.py* illustrate how the PomatoData can be used, which create data to model the German (DE) and central western Eurpoe (CWE) electricity systems.

The scope is defined by the settings arguments which is depicted below, here for the DE example. PomatoData will use power plant and network data for *grid_zones* and include all neighboring zones with all power plants, but network reduced to a single node. Weather year describes which year is used to derive the timeseries and capacity year defines which capacities to use. The time_horizon allows to create smaller datasets for a limited time horizon, e.g. May 2019. 

```python
    settings = {
        "grid_zones": ["DE"], 
        "weather_year": 2019,
        "capacity_year": 2030, 
        "co2_price": 60,
        "split_lines": True,
        "time_horizon": "01.05.2019 - 31.05.2019",
        }
```

The folder *pomato_datasets* contains an example dataset for Germany with the capacity year 2030.  

<img  src="https://github.com/richard-weinhold/pomato_data/blob/main/docs/_static/map_cwe.png?raw=true">


Installation
------------

To run PomatoData you should install it as a package. Clone the repository and run the *setup.py*
through *pip* in python. It is recommended to create a virtual environment and PomatoData pomato
into it, but not necessary. Execute the following commands from within the repository:
        
        python -m venv pomato_data
        ./pomato_data/Scripts/activate
        pip install -e .

Then you should be able to run the all the preprocessing, as described in *preprocess_input_data.py*
and subsequently create your own datasets, similar to the examples in *run_de_data.py* and
*run_cwe_data.py*.

Note, for *pip* on windows some packages can be difficult to install, as some dependencies are not
necessarily available. Those packages are specifically *Bottleneck*, *Rtree*, *rasterio*, *GDAL* and
*Fiona* which are dependencies for *geopandas* and *atlite*. Those can be installed individually
from the repository [Unofficial Windows Binaries for Python Extension Packages](https://www.lfd.uci.edu/~gohlke/pythonlibs/)
as *.whl* files. After the packages are installed, the PomatoData *setup.py* should successfully conclude.  


Release Status
--------------

PomatoData is part of the data pipeline i use in my PhD and is actively developed and will continue to change overtime. While there is a certain robustness to the data, especially since they represent great contributions by various parts of the open data community, the processing and combining is prone to errors both in program and concept, therefore I strongly recommend to use this data with caution and only after proper validation. 
