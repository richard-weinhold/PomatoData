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