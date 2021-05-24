# def fix_demand_file(wdir):
#     # Some files have hour in DateTime column
#     files = ["2020_1_ActualTotalLoad.csv"]
#     file_dir = wdir.joinpath("data_in/demand/entsoe_data")
    
#     for file in files:
#         demand_raw = pd.read_csv(file_dir.joinpath(file), header=0, encoding="UTF-16", sep="\t")
#         demand_raw['DateTime']  = pd.to_datetime(
#             demand_raw['Year'].astype(str) + '-' + \
#             demand_raw['Month'].astype(str) + '-' + \
#             demand_raw['Day'].astype(str) + ' ' + \
#             demand_raw['DateTime']
#             )

# def get_demand(wdir):
#     country_data, nuts_data = get_countries_regions_ffe()    
    
#     file_dir = wdir.joinpath("data_in/demand/entsoe_data")
#     files = [file for file in file_dir.glob("*.csv")]
    
#     usecols = ["DateTime", "ResolutionCode",
#                "AreaTypeCode", "AreaName", "MapCode", "TotalLoadValue"]
#     # Load Files    
#     demand = pd.DataFrame(columns=usecols)
#     for file in files:
#         ### Load Raw Data
#         demand_raw = pd.read_csv(file, header=0, encoding="UTF-16",
#                                sep="\t", usecols=usecols)
                
#         condition_country = (demand_raw.AreaTypeCode == "CTY") & \
#                             (demand_raw.MapCode.isin(country_data.index))
#         condition_60m = (demand_raw.ResolutionCode == "PT60M")
#         demand_raw = demand_raw.loc[condition_country&condition_60m]

#         demand = pd.concat([demand, demand_raw])
    
#     demand.sort_values("DateTime", inplace=True)
#     demand = demand.drop(["ResolutionCode", "AreaTypeCode"], axis=1)
#     # pcbf_df = pcbf_df.groupby(['DateTime', 'OutAreaName', 'OutMapCode', 'InAreaName', 'InMapCode']).sum().reset_index()
#     # pcbf_df = pcbf_df[~(pcbf_df.OutMapCode == pcbf_df.InMapCode)]
#     demand.columns = ["utc_timestep", "name", "zone", "demand"]
#     df = demand[demand.zone == "DK"]
#     return demand

def get_hydro_atlite(weather_year, cache_file_path, cache_file_name,
                     opsd_filepath, countries):

    weather_year = '2020'
    wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
    countries = ["NO"]

    cache_file_path = wdir.joinpath("data_temp")
    cache_file_name = "core"

        
    #Geoemtry shape for ERA5 cutout
    country_data, nuts_data = get_countries_regions_ffe()    
    
    x1, y1, x2, y2 = shapely.ops.cascaded_union(country_data.loc[countries, "geometry"].values).bounds
    #Path to store cutout
    # cutout_stor_path = 'data_temp\\'+cntr+'-'+weather_year
    cutout_stor_path = cache_file_path.joinpath(cache_file_name + '-' + str(weather_year))
           
    # Define cutout
    cutout = atlite.Cutout(path=str(cutout_stor_path),
                           module='era5',
                           x=slice(x1-.2, x2+.2), y=slice(y1-.2, y2+.2),
                           chunks={'time':100},
                           time=weather_year)
    cutout.prepare()
    
    cols = ["utc_timestamp", "nuts_id", "value"]
    wind = pd.DataFrame(columns=cols)
    pv = pd.DataFrame(columns=cols)
    for cntr in countries:
        zone = "BE"
        
        mato.data.plants["zone"] = mato.data.nodes.loc[mato.data.plants.node, "zone"].values
        tmp_plants = mato.data.plants[(mato.data.plants.technology == "ror")&(mato.data.plants.zone == zone)]
        
        # eps = 0.1
        # t = tmp_plants[(tmp_plants.lat < 40 )]
        
        # geometry = [shapely.geometry.Point(xy) for xy in zip(tmp_plants.lon, tmp_plants.lat)]
        # tmp_plants = gpd.GeoDataFrame(tmp_plants, crs="EPSG:4326", geometry=geometry)
        
        basins = gpd.read_file(r"C:\Users\riw\Documents\repositories\pomato_data\data_in\hydro\hybas_eu_lev12_v1c.zip")
        
        tmp = cutout.hydro(tmp_plants, basins, smooth=True)
        hydro_timeseries = tmp.to_pandas()
        hydro_timeseries = hydro_timeseries.T.fillna(0)
        hydro_timeseries.plot()
        
        # 1 m^3 = 1000 kg * 10 m 9.81 m/s^2 = 1000*10*9.81 Ws = (1000*10*9.81)/(3600 * 1e9) GWh
        
        hydro_timeseries_wh = hydro_timeseries*1000*10*9.81 / (3600*1e6)
        hydro_timeseries_wh.plot()
        
        # tmp_plants.loc[:, ["lat", "lon"]] =  53.060159, 8.865241 