
import sys 
import os
from pathlib import Path
import pandas as pd

if __name__ == "__main__":  

    # The purpose of this script is to preprocess the raw input data obtained 
    # from different sources and make them available to the PomatoData script, 
    # which combines the data. 

    # The idea is to preprocess data once (or few times), which enables to run 
    # PomatoData for different configuration, without starting from scratch. 
    
    # The folder structure in this repository is as follows:
    # - data_in: containing the raw data obtained from e.g. opsd or entso-e tp.
    # - data_tmp: contains temporary data or data obtained in the process of 
    #   preprocessing. This would be the atlite cutout file. 
    # - data_out: Contains the processed input data.

    # PomatoData only uses data from data out and filter/processes it into 
    # a single data-set compatible with POMATO. 
    
    if Path(os.path.abspath("")).name != "pomato_data":
        raise FileNotFoundError("Please Execute the script in the repository itself, use os.chdir() to change path")
    else: 
        wdir = Path(os.path.abspath(""))

    # %% Geographic information: FFE - Geodata 
    from pomato_data.auxiliary import get_countries_regions_ffe, get_eez_ffe
    
    # Onshore Areas are divided into NUTS3 areas
    zones, nuts_data = get_countries_regions_ffe(force_recalc=False)
    zones.to_csv(wdir.joinpath('data_out/zones/zones.csv'))
    nuts_data.to_csv(wdir.joinpath('data_out/zones/nuts_data.csv'))
    # Offshore areas are defined by the exclusive economic zones
    eez_region = get_eez_ffe(force_recalc=True)
    eez_region.drop("geometry", axis=1).to_csv(wdir.joinpath('data_out/zones/eez_wo_geometry.csv'))
    eez_region.to_csv(wdir.joinpath('data_out/zones/eez.csv'))
    
    # Note, both functions will cache previous downloads and return them unless
    # enforced by enforce_recalc=True
    
    # %% Weather Data: 
    # All weather data relies on atlite, it is therefore necessary to prepare 
    # a cutout including the desired region. In our case we want CWE plus 
    # neighboring countries. 

    # Process times are substantial. First, the cutout has to be made, with out 
    # geographic scope the final file takes 15GB or harddrive space, and 1-2 hours 
    # of download and process time. Once available and stored, re-instantiating the 
    # cutout takes seconds. 
    
    # Processing the timeseries from the cutout takes also significant amounts of 
    # time. With this scope it takes >60min
    
    from pomato_data.res import prepare_cutout, get_availabilities_atlite, offshore_eez_atlite

    opsd_filepath = wdir.joinpath("data_in/res/renewable_power_plants_DE.csv")
    countries = ["DE", "BE", "FR", "LU", "NL", "CH", "AT", "CZ", "DK", "PL", "SE", "ES", "PT", "UK", "NO", "IT"]
    weather_year = 2019
    cutout = prepare_cutout(weather_year, countries, wdir.joinpath("data_temp"), "core")
    # Wind Onshore, PV
    wind, pv = get_availabilities_atlite(cutout, countries, opsd_filepath)
    # Save Resulting Tables. 
    wind.to_csv(wdir.joinpath(f'data_out/res_availability/wind_availability_{weather_year}.csv'))
    pv.to_csv(wdir.joinpath(f'data_out/res_availability/pv_availability_{weather_year}.csv'))
    # Offshore
    offshore = offshore_eez_atlite(cutout)
    offshore.to_csv(wdir.joinpath(f'data_out/res_availability/offshore_availability_{weather_year}.csv'))
    
    # %% Hydro Plants with Inflows 
    # The inflows are calculated using atlite and the HydroSheds hydro basins 
    # data-set (https://www.hydrosheds.org/downloads) for EU, to determine inflows
    # we use levels 4-7, as indicated in the related publication. 
    # See: https://github.com/PyPSA/atlite/blob/1476c431360e05c0d51b4cc30248cb9e7294751e/atlite/convert.py#L529
    
    hydrobasins_path = wdir.joinpath("data_in/hydro/hydro_basins")
    from pomato_data.hydro import process_hydro_plants_with_atlite_inflows
    hydro_plants, inflows = process_hydro_plants_with_atlite_inflows(cutout, countries, hydrobasins_path)
    hydro_plants.to_csv(wdir.joinpath("data_out/hydro/plants.csv"))
    inflows.to_csv(wdir.joinpath(f"data_out/hydro/inflows_{weather_year}.csv"))
    
    # Storage Levels from ENTSO-E Transparency data. 
    
    from pomato_data.hydro import process_storage_level_entso_e
    storage_level = process_storage_level_entso_e(wdir, weather_year)
    storage_level.to_csv(wdir.joinpath(f"data_out/hydro/storage_level_{weather_year}.csv"))

    # %% RES Potentials
    # Potentials per NUTS3 area, from FFE open data portal 
    
    from pomato_data.res import get_potentials_ffe
    
    wind_potentials, pv_potentials = get_potentials_ffe()
    wind_potentials.to_csv(wdir.joinpath('data_out/res_potential/wind_potential.csv'))
    pv_potentials.to_csv(wdir.joinpath('data_out/res_potential/pv_potential.csv'))
    
    # %% RES capacities for wind and solar are determined based on weather year 
    # and zonal projected values from the FFE potentials. Other res without hydro 
    # is determined from OPSD renewable power plants data and mostly relates to biomass/biogas
    
    from pomato_data.res import other_res_capacities
    other_res_capacities = other_res_capacities(wdir)
    # Note the resulting table has NUTS 3 classification of year 2016, however 
    # all other data in this repo uses the FFE geo-data which is version 2013. 
    # The resulting file has to be converted using e.g. 
    # https://urban.jrc.ec.europa.eu/nutsconverter/#/
    other_res_capacities.to_csv(wdir.joinpath('data_out/res_capacity/other_res_v16.csv'))

    # The resulting file is included with this repository. 
    
    # %% Network (Nodes&Lines (lines contain AC/DC and transformer branches))
    # This requires the GridKit data from pyPSA downloaded into data_in/GridKit
    # This currently uses the version from jan 2020. 
    
    from pomato_data.grid import process_gridkit_data
    gridkit_filepath = wdir.joinpath("data_in/GridKit")
    nodes, lines = process_gridkit_data(gridkit_filepath)
        
    # Manually adding additionally dc-lines 
    add_dclines = pd.read_csv(wdir.joinpath("data_in/grid/add_dclines.csv"), index_col=0)
    tmp_lines = pd.concat([lines, add_dclines], axis=0)
    tmp_lines.to_csv(wdir.joinpath("data_out/lines/lines.csv"))
    nodes.to_csv(wdir.joinpath("data_out/nodes/nodes.csv"))
    
    # %% Exchange: commercial's and physical 
    # For commercial exchange, this requires the ScheduledCommercialExchanges
    # for the chosen year in data_in\exchange\commercial_exchange and for physical 
    # cross-boarder flow CrossBorderPhysicalFlow. 
    # Both are available from the ENTSO-E Transparency platform FTP server. 
    from pomato_data.exchange import (process_commercial_exchange_entso_e, 
                                      process_physical_crossborder_flow_entso_e)
                                     
    physical_crossborder_flow = process_physical_crossborder_flow_entso_e(wdir, weather_year)
    commercial_exchange = process_commercial_exchange_entso_e(wdir, weather_year)
    physical_crossborder_flow.to_csv(wdir.joinpath(f"data_out/exchange/physical_crossborder_flow_{weather_year}.csv"))
    commercial_exchange.to_csv(wdir.joinpath(f"data_out/exchange/commercial_exchange_{weather_year}.csv"))

    
    # %% Demand:
    # Derived from ActualTotalLoad available from the ENTSO-E TP FTP server. 
    # The monthly csv files have to be available in data_in\demand\ensto-e
    from pomato_data.demand import get_demand_entso_e
    demand = get_demand_entso_e(wdir, weather_year)
    demand.to_csv(wdir.joinpath(f'data_out/demand/demand_{weather_year}.csv'))
    
    # To distribute zonal demand onto nodes, we use standard load profiles from BDEW 
    # and nuts3 regionalized GDP data from Eurostat, both are freely available 
    # but packages with the repo to ensure compatibility. 
    # Extract SLP and GDP data if not already in data_in/demand
    gdp_slp_folder = wdir.joinpath("data_in/demand/gdp_data")
    if not gdp_slp_folder.is_dir():
        import zipfile
        with zipfile.ZipFile(gdp_slp_folder.with_suffix(".zip"), 'r') as zip_ref:
            zip_ref.extractall(gdp_slp_folder.parent)
    
    # %% Conventional Power Plants 
    
    # Conventional power plants come from the OPSD package, but is processed to 
    # include coordinates for (almost) all plants. This processed file is also 
    # packaged with this repo and then further processed to include efficuency 
    # estimates. 
    from pomato_data.plants import process_plants 
    plants = process_plants(wdir)
    plants.to_csv(wdir.joinpath('data_out/plants/plants.csv'))
    
    
    
    
  