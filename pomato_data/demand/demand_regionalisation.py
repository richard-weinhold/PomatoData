import os
import pyproj

import geopandas as gpd
from shapely.geometry import Point
import pandas as pd

def nodal_demand(ddir, countries, demand_el, nodes):
    nodal_demand = pd.DataFrame()
    for country in countries:
        print(country)
        nodal_demand = pd.concat([nodal_demand, get_nodal_demand(ddir, country, demand_el, nodes)], axis=1)

    return nodal_demand


def only_nuts3_zones(ddir):
    nuts3 = pd.read_csv(ddir.joinpath('data_in/demand/gdp_data/NUTS_AT_2013.csv'),
                    encoding = "ISO-8859-1", index_col=0,
                    usecols=["CNTR_CODE", "NUTS_ID", "NUTS_NAME"])
    ### Remove Wierd things and NANs
    nuts3.NUTS_ID = nuts3.NUTS_ID.str.replace("-", "")
    nuts3.NUTS_ID = nuts3.NUTS_ID.str.upper()

    nuts3 = nuts3[~nuts3.NUTS_ID.isnull()]
    nuts3 = nuts3[nuts3.NUTS_ID.apply(lambda x: (len(x) == 5) and (sum(c.isdigit() for c in x) in [0,1,2,3]))]
    return nuts3

def standard_load_profiles(ddir):
    # ## Read BDEW Standard Load Profiles for household and Commercial
    # Source: https://www.stromnetz.berlin/de/netznutzer.htm
    slp_h0_60min = pd.read_excel(ddir.joinpath('data_in/demand/gdp_data/bdew-slp-H0-2015.xls'),
                             sheet_name='Haushaltsprofil 2015', usecols=[0,2],
                             index_col=0, names=["index", 'load']).resample('1h').sum()
    slp_h0_60min.index = slp_h0_60min.index.round(freq='60min')

    slp_g0_60min = pd.read_excel(ddir.joinpath('data_in/demand/gdp_data/bdew-slp-G0-2015.xls'),
                             sheet_name='Gewerbe allgemein 2015', usecols=[0,2],
                             index_col=0, names=["index", 'load']).resample('1h').sum()
    slp_g0_60min.index = slp_g0_60min.index.round(freq='60min')

    return slp_h0_60min, slp_g0_60min

def read_comsumption_data(ddir):
    
    consumtion = pd.read_csv(ddir.joinpath('data_in/demand/gdp_data/nrg_105a_1_Data.csv'),
                             encoding = "ISO-8859-1", usecols=["GEO",'INDIC_NRG', 'Value'])
    consumtion.loc[:, "Value"] = pd.to_numeric(consumtion.Value.str.replace(".","", regex=False), errors='coerce').fillna(0).values
    
    consumtion = consumtion.pivot(index="GEO", columns="INDIC_NRG", values="Value")
    consumtion.columns = ["commercial", "household", "misc", "misc_other", "industry", "total"]

    return consumtion

def load_nuts_data(ddir):

    # Population from EUROSTAT
    # Source: http://ec.europa.eu/eurostat/de/web/rural-development/data
    pop_nuts3 = pd.read_csv(ddir.joinpath('data_in/demand/gdp_data/nama_10r_3popgdp_1_Data.csv'),
                        encoding = "ISO-8859-1", usecols=["GEO", 'Value'], index_col=0)

    pop_nuts3.Value = pd.to_numeric(pop_nuts3.Value.str.replace(",",""), errors='coerce').fillna(0)
    pop_nuts3.Value *= 1000
    pop_nuts3.columns = ["pop"]

    gva_nuts3 = pd.read_csv(ddir.joinpath('data_in/demand/gdp_data/nama_10r_3gva_1_Data.csv'),
                        usecols=["GEO",'NACE_R2', 'Value'], index_col=0)
    # Gross Values Added
    gva_nuts3.Value = pd.to_numeric(gva_nuts3.Value.str.replace(",",""), errors='coerce').fillna(0)
    gva_nuts3.Value[gva_nuts3.NACE_R2 == "B-E"] + gva_nuts3.Value[gva_nuts3.NACE_R2 == "F"]
    gva_nuts3 = gva_nuts3.pivot(columns="NACE_R2", values="Value")
    gva_nuts3["gva_industry"] = gva_nuts3["B-E"] + gva_nuts3["F"]
    gva_nuts3["gva_other"] =  gva_nuts3["TOTAL"] - gva_nuts3["gva_industry"]
    gva_nuts3 = gva_nuts3[["gva_other", "gva_industry"]]

    return pop_nuts3, gva_nuts3

def get_nodal_demand(ddir, country, demand, nodes, scaling=None):
    
    # %%
    # country = "DK"
    # scaling = None
    # ddir = self.wdir
    # nodes = self.nodes
    # demand = self.demand

    # %%    
    print("Regionalizing demand for", country)

    node_data = nodes[(nodes.zone == country) & (nodes.demand)]
    node_data.index.name = "id"
    node_data.reset_index(inplace=True)
    geometry = [Point(xy) for xy in zip(node_data.lon, node_data.lat)]
    node_data = gpd.GeoDataFrame(node_data, crs="EPSG:4326", geometry=geometry)
    

    ## 2. Using ENTSO-E Load Time Series from OPSD
    demand_el = demand[country].to_frame()
    demand_el.columns = ["load"]
    
    if len(node_data) == 1:
        demand_el.columns = [node_data.id[0]]
        return demand_el

    nuts3 = only_nuts3_zones(ddir)
    pop_nuts3, gva_nuts3 = load_nuts_data(ddir)
    consumtion = read_comsumption_data(ddir)

    slp_h0_60min, slp_g0_60min = standard_load_profiles(ddir)

    slp_h0_60min.index = slp_h0_60min.index.map(lambda t: t.replace(year=demand_el.index[0].year))
    slp_h0_60min = slp_h0_60min[slp_h0_60min.index.isin(demand_el.index)]
    slp_g0_60min.index = slp_g0_60min.index.map(lambda t: t.replace(year=demand_el.index[0].year))
    slp_g0_60min = slp_g0_60min[slp_g0_60min.index.isin(demand_el.index)]

    split_household = consumtion[consumtion.index == country].household.values[0]
    split_commercial = consumtion[consumtion.index == country].commercial.values[0]
    # split_industry = consumtion[consumtion.index == country].industry.values[0]
    split_total = consumtion[consumtion.index == country].total.values[0]

    if isinstance(scaling, dict):
        hh = f"Household {round(scaling['household']*100, 2)}%"
        com = f"Commercial {round(scaling['commercial']*100, 2)}%"
        ind = f"Industry {round(scaling['industry']*100, 2)}%"
        print("Splits between Sectors from Input: Household " + hh + ", " + com + ", " + ind)
        scale_household = scaling['household']*(demand_el.load.sum()/(slp_h0_60min.load.sum()))
        scale_commercial = scaling['commercial']*(demand_el.load.sum()/(slp_g0_60min.load.sum()))
        print(f"Scaling demand by {scaling['demand']*100}%")
        demand_el = demand_el*scaling['demand']
    else:
        hh = f"Household {round(split_household/split_total*100, 2)}%"
        com = f"Commercial {round(split_commercial/split_total*100, 2)}%"
        ind = f"Industry {round((1 - split_household/split_total - split_commercial/split_total)*100, 2)}%"
        print("Splits between Sectors from Input: Household " + hh + ", " + com + ", " + ind)
        scale_household = (split_household/split_total)*(demand_el.load.sum()/(slp_h0_60min.load.sum()))
        scale_commercial = (split_commercial/split_total)*(demand_el.load.sum()/(slp_g0_60min.load.sum()))

    # ## 3.1 Scale to AGEB Energy Balance
    loadprofile_household = slp_h0_60min  * scale_household
    loadprofile_commercial = slp_g0_60min * scale_commercial

    demand_el.index = demand_el.index.astype('datetime64[ns]')
    loadprofile_industry = demand_el - loadprofile_household - loadprofile_commercial
    loadprofile_industry[loadprofile_industry['load'] < 0] = 0
    # loadprofile_industry.sum()

    pop_nuts3 = pop_nuts3[pop_nuts3.index.isin(nuts3.NUTS_ID[nuts3.index == country])]
    gva_nuts3 = gva_nuts3[gva_nuts3.index.isin(nuts3.NUTS_ID[nuts3.index == country])]

    ## Merge all NUTS3 stats
    stats_nuts3 = pd.merge(gva_nuts3,pop_nuts3,left_index=True,right_index=True)

    ## 4.4 Calculate Scale Factors
    sf_nuts3 = pd.concat([stats_nuts3['pop']/stats_nuts3['pop'].sum(),
                          stats_nuts3['gva_industry']/stats_nuts3['gva_industry'].sum(),
                          stats_nuts3['gva_other']/stats_nuts3['gva_other'].sum()],axis=1)

    # sf_nuts3.head()

    #Regionalize Load Profiles to NUTS3
    lp_household_reg = pd.DataFrame([sf_nuts3['pop'] for n in range(len(loadprofile_household['load'].values))],index=loadprofile_household.index)
    lp_household_reg = lp_household_reg.mul(loadprofile_household['load'],axis=0)

    lp_commercial_reg = pd.DataFrame([sf_nuts3['gva_other'] for n in range(len(loadprofile_commercial['load'].values))],index=loadprofile_commercial.index)
    lp_commercial_reg = lp_commercial_reg.mul(loadprofile_commercial['load'],axis=0)

    lp_industry_reg = pd.DataFrame([sf_nuts3['gva_industry'] for n in range(len(loadprofile_industry['load'].values))],index=loadprofile_industry.index)
    lp_industry_reg = lp_industry_reg.mul(loadprofile_industry['load'],axis=0)

    # lp_reg = pd.concat({'household':lp_household_reg,
    #                     'commercial':lp_commercial_reg,
    #                     'industry':lp_industry_reg},
    #                       axis=1)

    # NUTS3 Regions
    # Source: http://ec.europa.eu/eurostat/de/web/gisco/geodata/reference-data/administrative-units-statistical-units/nuts
    nuts_rg_01m_2013 = gpd.read_file(str(ddir.joinpath('data_in/demand/gdp_data/NUTS_2013_01M_SH/data/NUTS_RG_01M_2013.shp')))
    nuts_rg_01m_2013 = nuts_rg_01m_2013.loc[(nuts_rg_01m_2013['STAT_LEVL_'] == 3) & (nuts_rg_01m_2013['NUTS_ID'].str.contains(country, regex=True))]
    # nuts_rg_01m_2013_de.set_index('NUTS_ID',inplace=True)
    nuts_rg_01m_2013 = nuts_rg_01m_2013.to_crs('epsg:4326')


    # ## 6.3 Plot NUTS3 and Nodes
    # base = nuts_rg_01m_2013.plot(color='white',figsize=(15, 20))
    # node_data.plot(ax=base, marker='o', color='red', markersize=5)

    # %%
    # # 7. Distribute Load to Nodes
    # ## 7.1 Spacial Joins
    # Spacial join of NUTS regions to each node.
    nuts_to_nodes = gpd.sjoin(node_data, nuts_rg_01m_2013, how='left', op='within')
    
    # node_data.iloc[0].geometry.within(nuts_rg_01m_2013.iloc[0].geometry)
    
    if not nuts_to_nodes[nuts_to_nodes.NUTS_ID.isnull()].empty:
        print(nuts_to_nodes.loc[nuts_to_nodes.NUTS_ID.isnull(), ["id", "lat", "lon"]])
        print("WARNING: Some Nodes with demand boolean are not in NUTS Areas, \
               possibly outside of coutry or in sea area")
    
    nuts_to_nodes["weighting"] = 0.5
    nuts_to_nodes["weighting"] = nuts_to_nodes["voltage"].map({380: 2, 220: 1, 132: 0.5})  
    if any(nuts_to_nodes.weighting.isna()):
        nuts_to_nodes["weighting"] = 1
    # nuts_to_nodes["weighting"] = nuts_to_nodes["voltage"].map({380: 1, 220: 1, 132: 1})  

    # Spacial join of nodes to each NUTS region.
    node_in_nuts = gpd.sjoin(nuts_rg_01m_2013, node_data, how='left', op='contains')

    # ## 7.2 Regions Containing (Multiple) Nodes
    # Count the number of nodes within one region.
    nuts_multiple_nodes = nuts_to_nodes[["NUTS_ID", "weighting"]].groupby('NUTS_ID').count()
    nuts_multiple_nodes["sum_weight"] = nuts_to_nodes[["NUTS_ID", "weighting"]].groupby('NUTS_ID').sum()
    nuts_to_nodes["weight"] = (nuts_to_nodes.loc[:, "weighting"]/nuts_multiple_nodes.loc[nuts_to_nodes['NUTS_ID'], "sum_weight"].values)

    # Distribute the load equally to all nodes in one region.
    lp_household_node = lp_household_reg[nuts_to_nodes['NUTS_ID'].values].multiply(nuts_to_nodes["weight"].values)
    lp_household_node.columns = nuts_to_nodes.index

    lp_commercial_node = lp_commercial_reg[nuts_to_nodes['NUTS_ID'].values].multiply(nuts_to_nodes["weight"].values)
    lp_commercial_node.columns = nuts_to_nodes.index

    lp_industry_node = lp_industry_reg[nuts_to_nodes['NUTS_ID'].values].multiply(nuts_to_nodes["weight"].values)
    lp_industry_node.columns = nuts_to_nodes.index

    # %%
    # ## 7.3 Add Regions not Containing any Nodes
    # Find NUTS regions without any nodes within and calculate their centroids.
    nuts_no_node = node_in_nuts[node_in_nuts["index_right"].isna()]
    # nuts_centroids = nuts_rg_01m_2013.to_crs('epsg:2953').centroid.to_crs("epsg:4326")
    # nuts_no_node = nuts_centroids.loc[nuts_no_node.index].to_crs('epsg:2953')    
    nuts_centroids = nuts_rg_01m_2013.centroid
    nuts_no_node = nuts_centroids.loc[nuts_no_node.index]

    # Find the closest node to the region's centroid and map them
    dist = []
    for j in range(len(node_data)):
        dist.append(nuts_no_node.distance(node_data.geometry.iloc[j]))
    dist = pd.DataFrame(dist)

    nuts_no_node_map = pd.DataFrame(dist.idxmin(),columns=['node_index'])
    nuts_no_node_map = nuts_no_node_map.join(pd.DataFrame(nuts_rg_01m_2013['NUTS_ID']))

    nuts_no_node_map = nuts_no_node_map.reset_index()
    nuts_no_node_map.rename(columns={'index': 'nuts_index'}, inplace=True)

    # Add those regions to the nodal load patterns
    for i in nuts_no_node_map.index:
        lp_household_node.loc[:, nuts_no_node_map['node_index'].loc[i]] += lp_household_reg.loc[:, nuts_no_node_map['NUTS_ID'].loc[i]]
        lp_commercial_node.loc[:, nuts_no_node_map['node_index'].loc[i]] += lp_commercial_reg[nuts_no_node_map['NUTS_ID'].loc[i]]
        lp_industry_node.loc[:, nuts_no_node_map['node_index'].loc[i]] += lp_industry_reg[nuts_no_node_map['NUTS_ID'].loc[i]]

    lp_node = (lp_household_node + lp_commercial_node + lp_industry_node)
    lp_node.columns = node_data.id

    # Check whether the total load still corresponds with AGEB energy balance
    print("Calc", lp_node.sum().sum())
    print("TS", demand_el.load.sum())

    # Barchart
    # barplot = pd.DataFrame(index=["Demand"], columns=["Household", "Commercial", "Industry"])
    # barplot.loc["Demand", "Household"] = lp_household_node.sum().sum()/1e6
    # barplot.loc["Demand", "Commercial"] = lp_commercial_node.sum().sum()/1e6
    # barplot.loc["Demand", "Industry"] = lp_industry_node.sum().sum()/1e6
    # barplot.plot.bar(stacked=True)

    # Plot the ENTSO-E profile and the aggregated load profile
    # base = demand_el.rename(columns={'load': 'ENTSO-E'}).plot(color='blue',figsize=(20, 10))
    # pd.DataFrame(lp_node.sum(axis=1)).rename(columns={0: 'Load Profile Total'}).plot(ax = base, color='red',figsize=(20, 10), legend = 'Load Profile Total');
    # pd.DataFrame(lp_household_node.sum(axis=1)).rename(columns={0: 'Load Profile Household'}).plot(ax = base, color='green',figsize=(20, 10), legend = 'Load Profile Household');
    # pd.DataFrame(lp_commercial_node.sum(axis=1)).rename(columns={0: 'Load Profile Commercial'}).plot(ax = base, color='pink',figsize=(20, 10), legend = 'Load Profile Commercial');
    # pd.DataFrame(lp_industry_node.sum(axis=1)).rename(columns={0: 'Load Profile Industry'}).plot(ax = base, color='purple',figsize=(20, 10), legend = 'Load Profile Industry');
    # %%
    return lp_node






