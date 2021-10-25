

import pandas as pd
from pathlib import Path
from pomato_data.auxiliary import get_countries_regions_ffe
   
def get_demand_entso_e(wdir, year):
    country_data, nuts_data = get_countries_regions_ffe()    
    usecols = ["DateTime", "AreaTypeCode", "MapCode", "TotalLoadValue"]
    
    file_dir = wdir.joinpath("data_in/demand/ensto-e")

    files = [file for file in file_dir.glob("*.csv") if str(year) in str(file)]

    # Load Files    
    demand = pd.DataFrame()
    for file in files:
        ### Load Raw Data
        try: 
            demand_raw = pd.read_csv(file, header=0, encoding="UTF-16",
                                     sep="\t", usecols=usecols)
        except UnicodeError:
            demand_raw = pd.read_csv(file, header=0, sep="\t", usecols=usecols)

        demand_raw["MapCode"] = demand_raw["MapCode"].replace("GB", "UK")
        demand_raw.DateTime = pd.to_datetime(demand_raw.DateTime).astype('datetime64[ns]')
        demand_raw['DateTime'] = demand_raw['DateTime'].apply(lambda x:x.replace(minute=0))
        demand_raw = demand_raw.groupby(usecols[:-1]).mean().reset_index()
        
        condition_cty = (demand_raw.AreaTypeCode == "CTY") 
        condition_cta = (demand_raw.AreaTypeCode == "CTA")
        condition_2 = demand_raw.MapCode.isin(country_data.index)
                    
        demand_cty = demand_raw.loc[condition_cty & condition_2].copy()
        demand_cta = demand_raw.loc[condition_cta & condition_2].copy()
        
        cols = ["DateTime", "MapCode", "TotalLoadValue"]
        demand_raw = pd.merge(demand_cty[cols], demand_cta[cols], on=cols[:-1], how="outer", suffixes=("", "_cta"))
        condition = (demand_raw.TotalLoadValue.isna())&(demand_raw.TotalLoadValue_cta.notna())

        demand_raw.loc[condition, "TotalLoadValue"] = demand_raw.loc[condition, "TotalLoadValue_cta"]     
        demand_raw = demand_raw.drop("TotalLoadValue_cta", axis=1)
        demand = pd.concat([demand, demand_raw])
        
    demand.sort_values("DateTime", inplace=True)

    # Fix missing values with the one an hour before
    demand = demand.pivot(index="DateTime", columns = "MapCode", values="TotalLoadValue")
    for zone in demand.columns:
        for i in demand.loc[demand[zone].isna(), zone].index:
               demand.iloc[demand.index.get_loc(i)][zone] = demand.iloc[demand.index.get_loc(i) - 1][zone]
    demand = demand.reset_index().rename(columns={"DateTime": "utc_timestamp"})
    
    return demand 

if __name__ == "__main__": 
    import pomato_data

    wdir = Path(pomato_data.__path__[0]).parent 
    year = 2015
    demand = get_demand_entso_e(wdir, year)
    demand.to_csv(wdir.joinpath(f'data_out/demand/demand_{year}.csv'))
    
