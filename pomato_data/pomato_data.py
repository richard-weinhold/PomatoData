


import os

import requests
import pandas as pd
from pathlib import Path
import shapely
from shapely.geometry import Point, LineString
import pyproj
import geopandas as gpd
import numpy as np 
import datetime as dt
import itertools  
from scipy import sparse
import shutil

os.chdir(r'C:\Users\riw\Documents\repositories\pomato_data')
from pomato_data.auxiliary import get_countries_regions_ffe, distance, \
    load_data_structure, add_timesteps, match_plants_nodes
from pomato_data.res.capacities import regionalize_res_capacities, existing_offshore_wind_capacities
from pomato_data.demand.demand_regionalisation import nodal_demand


# %%
class PomatoData():
    
    def __init__(self, wdir, settings):
        self.wdir = wdir
        self.settings = settings
        self.data_structure = load_data_structure(self.wdir)

        startend = self.settings["time_horizon"].split(" - ")
        self.time_horizon = {"start": dt.datetime.strptime(startend[0], "%d.%m.%Y"),
                             "end": dt.datetime.strptime(startend[1], "%d.%m.%Y")}
        
        self.load_data()
        self.process_lines_nodes()
        self.split_double_lines()
        self.connect_small_subnetworks()

        self.process_zones()
        self.process_demand()
        self.process_plants()
        self.process_res_plants()
        self.marginal_costs()
        self.fix_plant_connections()

        self.process_availabilites()
        self.process_offshore_plants()
        self.create_basic_ntcs()

    def load_data(self):
        
        self.zones = pd.read_csv(self.wdir.joinpath("data_out/zones/zones.csv"), index_col=0)
        self.nuts_data = pd.read_csv(self.wdir.joinpath("data_out/zones/nuts_data.csv"), index_col=0)
        self.nodes = pd.read_csv(self.wdir.joinpath("data_out/nodes/nodes.csv"), index_col=0)
        self.lines = pd.read_csv(self.wdir.joinpath("data_out/lines/lines.csv"), index_col=0)

        self.plants = pd.read_csv(self.wdir.joinpath("data_out/plants/plants.csv"), index_col=0)
        self.fuel = pd.read_csv(self.wdir.joinpath("data_out/fuel/fuel.csv"), index_col=0)
        self.technology = pd.read_csv(self.wdir.joinpath("data_out/technology/technology.csv"), index_col=0)
        
        self.demand_el = pd.read_csv(self.wdir.joinpath('data_out/demand/demand.csv'), index_col=0)
        self.demand_el.utc_timestamp = pd.to_datetime(self.demand_el.utc_timestamp).astype('datetime64[ns]')
        
    def process_zones(self):
        # Reset zones
        self.zones = self.zones[self.zones.index.isin(self.nodes.zone)]
        non_grid_zones = [z for z in self.zones.index if z not in self.settings["grid_zones"]]
        add_nodes = []
        add_dclines = []
        geod = pyproj.Geod(ellps="WGS84")
        for z in non_grid_zones:
            # z = "NO"
            if z == "NO":
                lat, lon = 63.342806, 10.459677
            else:
                lon, lat = shapely.wkt.loads(self.zones.loc[z, "geometry"]).centroid.coords[0]
            name = self.zones.loc[z, "name"]
            # substation', 'voltage', 'name', 'lat', 'lon', 'zone', 'info', 'demand', 'slack'
            add_nodes.append(["n" + z, "n" + z, 500, name, lat, lon, z, "", True, True])
            self.nodes.loc[self.nodes.zone == z, "demand"] = False
            # 'node_i', 'node_j', 'name_i', 'name_j', 'voltage', 'length', 'under_construction', 
            # 'geometry', 'technology', 'capacity', 'status', 'commissioned'
            for n in self.nodes[self.nodes.zone == z].index:
                geometry = LineString([Point(self.nodes.loc[n, ["lon", "lat"]]), Point(lon, lat)])
                length = geod.geometry_length(geometry)
                add_dclines.append(["dc" + z + "_" + n, "n" + z, n, z, n, 500, length, False,
                                    geometry.wkt, "dc", 1e5, "online", 1900])
        
        add_nodes = pd.DataFrame(index=np.array(add_nodes)[:, 0], columns=self.nodes.columns, data=np.array(add_nodes)[:, 1:])
        add_dclines = pd.DataFrame(index=np.array(add_dclines)[:, 0], columns=self.dclines.columns, data=np.array(add_dclines)[:, 1:])
        
        self.nodes = self.nodes.append(add_nodes)
        self.dclines = self.dclines.append(add_dclines)

        self.dclines[["capacity", "length"]] = self.dclines[["capacity", "length"]].astype("float")
        self.nodes[["lat", "lon", "voltage"]] = self.nodes[["lat", "lon", "voltage"]].astype("float")
        # self.nodes.loc[:, ["lat", "lon"]].dtypes
        
    def process_lines_nodes(self):
        """Process Nodes and Lines Data."""
        nodes_in_co = self.nodes.index[self.nodes.zone.isin(self.settings["grid_zones"])]
        self.lines = self.lines[(self.lines.node_i.isin(nodes_in_co)) |
                                (self.lines.node_j.isin(nodes_in_co))]
        # Filter Online
        # condition_online = (self.lines.status == "online") | (self.lines.commissioned <= self.settings["year"])
        condition_online = (self.lines.node_i.isin(self.nodes.index))
        # 
        self.dclines = self.lines[(self.lines.technology == "dc") & condition_online]
        self.lines = self.lines[(self.lines.technology.isin(["ac_line", "ac_cable", "transformer"]) & condition_online)]
        # Only use nodes with lines attached
        self.nodes = self.nodes[(self.nodes.index.isin(self.lines.node_i)) |
                                (self.nodes.index.isin(self.lines.node_j)) |
                                (self.nodes.index.isin(self.dclines.node_i)) |
                                (self.nodes.index.isin(self.dclines.node_j))]

        self.nodes["slack"] = False
        self.lines["contingency"] = True
        self.lines.loc[self.lines.technology == "transformer", "contingency"] = False
        self.dclines = self.dclines.drop(["x", "r", "x_pu", "r_pu", "circuits"], axis=1)
        

    def split_double_lines(self):
        """Split Lines with circuits > 1 into multiples."""
        self.lines = self.lines.reset_index()
        self.lines["part of # circuits"] = self.lines.circuits
        for nr_c in range(2, int(self.lines.circuits.max()) + 1):
            for i in range(2, nr_c + 1):
                tmp_df = self.lines[self.lines.circuits == nr_c].copy()
                tmp_df.loc[:, "index"] = tmp_df.loc[:, "index"] + "_" + str(i) 
                tmp_df.loc[:, "circuits"] = 1
                self.lines = self.lines.append(tmp_df)
            self.lines.loc[self.lines.circuits == nr_c, "index"] = self.lines.loc[self.lines.circuits == nr_c, "index"] + "_1"
            self.lines.loc[self.lines.circuits == nr_c, "circuits"] = 1
        self.lines = self.lines.set_index("index")
        self.lines.drop(["circuits"], axis=1, inplace=True)
    
    
    def fix_plant_connections(self):
        """Make sure nodes are connected with enough line capacity"""
        tmp_nodes = self.nodes.copy()
        for node in tmp_nodes.index:
            self.dclines.loc[(node == self.dclines.node_i)|(node == self.dclines.node_j), "capacity"].sum() 
            self.dclines.loc[(node == self.dclines.node_i)|(node == self.dclines.node_j), "capacity"]
            
        line_node_capacity = [self.dclines.loc[(node == self.dclines.node_i)|(node == self.dclines.node_j), "capacity"].sum() 
                              + self.lines.loc[(node == self.lines.node_i)|(node == self.lines.node_j), "capacity"].sum() for node in tmp_nodes.index]

        tmp_nodes["capacity"] = self.plants[["g_max", "node"]].groupby("node").sum()
        tmp_nodes["capacity"] = tmp_nodes["capacity"].fillna(0)
        tmp_nodes["demand_max"] = self.demand_el.drop("utc_timestamp", axis=1).max(axis=0)
        tmp_nodes["demand_min"] = self.demand_el.drop("utc_timestamp", axis=1).min(axis=0)
        tmp_nodes["line_capacity"] = line_node_capacity
        tmp_nodes["needed_capacity_min"] = tmp_nodes["demand_max"] - tmp_nodes["capacity"]
        tmp_nodes["needed_capacity_max"] = tmp_nodes["capacity"] - tmp_nodes["demand_min"]
        
        tmp_nodes["needed_capacity"] =tmp_nodes[["needed_capacity_min", "needed_capacity_max"]].max(axis=1)
        tmp_nodes = tmp_nodes[tmp_nodes.needed_capacity > tmp_nodes["line_capacity"]]        
        for node in tmp_nodes.index:
            needed_capacity = tmp_nodes.loc[node, "needed_capacity"]
            line_cap = tmp_nodes.loc[node, "line_capacity"]
            print(f"Node {node} needs {needed_capacity}MW of transport capacity but only {line_cap}MW of line capacity connected. Increasing capacity accordingly.")
            self.lines.loc[(node == self.lines.node_i)|(node == self.lines.node_j), "capacity"] *= needed_capacity/line_cap*1.001

    def process_plants(self):
        
        self.plants = self.plants[self.plants.zone.isin(self.zones.index)].reset_index(drop=True)
        self.plants = self.plants[~self.plants.status.isin(["shutdown", "shutdown_temporary"])]
        self.plants["heatarea"] = None
        self.plants.index = "p" + self.plants.index.astype(str)
        self.plants.loc[self.plants.g_max.isnull(), "g_max"] = 0
        self.plants.loc[self.plants.h_max.isnull(), "h_max"] = 0
        
        self.plants.loc[self.plants.plant_type.isin(["hydro_res", "hydro_psp"]),
                        "storage_capacity"] = self.plants.g_max * 24*2
        
        self.plants = self.plants[(self.plants.g_max>0)]
        self.plants["node"] = None
        self.plants = match_plants_nodes(self.plants, self.nodes)
        
    def marginal_costs(self):
        """Calculate the marginal costs for plants that don't have it manually set.

        Marginal costs are calculated by: mc = fuel price / eta + o&m + CO2 costs
        """
        # self = data
        co2_price = self.settings["co2_price"]
        tmp_costs = pd.merge(self.fuel, self.technology[['fuel', 'technology', "plant_type", 'variable_om']],
                             how='left', left_index=True, right_on=['fuel'])
        self.plants[["mc_el", "mc_heat"]] = np.nan
        tmp_plants = self.plants.copy()
        
        tmp_plants = pd.merge(tmp_plants, tmp_costs, how='left', on=['technology', 'fuel', "plant_type"])
        tmp_plants.mc_el = tmp_plants.fuel_price / tmp_plants.eta + tmp_plants.variable_om + \
                           tmp_plants.co2_content * co2_price
            
        # t = tmp_plants.loc[tmp_plants.mc_el.isna(), ['fuel', 'technology', "plant_type", 'variable_om', "fuel_price", "co2_content", "eta", "mc_el"]]
        self.plants.loc[:, "mc_el"] = tmp_plants.mc_el.values
        self.plants.loc[:, "mc_heat"] = 0
        
    def process_demand(self):

        self.demand_el = self.demand_el[(self.demand_el.utc_timestamp >= self.time_horizon["start"]) & \
                                  (self.demand_el.utc_timestamp < self.time_horizon["end"])]
        self.demand_el = self.demand_el.set_index("utc_timestamp")
        
        self.demand_el = nodal_demand(self.wdir, list(self.zones.index), self.demand_el, self.nodes)
        # Set demand of missing Nodes to 0
        for node in [n for n in self.nodes.index if n not in self.demand_el.columns]:
            self.demand_el[node] = 0
        self.demand_el.index.name = "utc_timestamp"
        self.demand_el = add_timesteps(self.demand_el)
    
    def process_res_plants(self):
        res_plants = regionalize_res_capacities(self.wdir, self.nodes.copy(), self.zones.index, self.technology)
        self.plants = pd.concat([self.plants, res_plants])
        self.plants.loc[self.plants.plant_type.isin(["hydro_res", "hydro_psp"]),
                "storage_capacity"] = self.plants.g_max * 24*2
        
    def process_availabilites(self):
        
        plants = self.plants[self.plants.technology.isin(["solar", "wind onshore"])]
        nodes =  self.nodes[self.nodes.index.isin(plants.node)].copy()
        
        wind_availability = pd.read_csv(self.wdir.joinpath('data_out/res_availability/wind_availability.csv'), index_col=0)
        pv_availability = pd.read_csv(self.wdir.joinpath('data_out/res_availability/pv_availability.csv'), index_col=0)

        wind_availability = wind_availability.pivot(index="utc_timestamp", columns="nuts_id", values="value")
        wind_availability.index = pd.to_datetime(wind_availability.index).astype('datetime64[ns]')
        wind_availability = wind_availability.sort_index()
        
        pv_availability = pv_availability.pivot(index="utc_timestamp", columns="nuts_id", values="value")
        pv_availability.index = pd.to_datetime(pv_availability.index).astype('datetime64[ns]')
        pv_availability = pv_availability.sort_index()
        
        wind_availability = wind_availability[(wind_availability.index >= self.time_horizon["start"]) & \
                                              (wind_availability.index < self.time_horizon["end"])]
        
        pv_availability = pv_availability[(pv_availability.index >= self.time_horizon["start"]) & \
                                          (pv_availability.index < self.time_horizon["end"])]
        
        country_data, nuts_data = get_countries_regions_ffe()    
        # nuts_data['geometry'] = nuts_data['geometry'].apply(shapely.wkt.loads)
        nuts_data = gpd.GeoDataFrame(nuts_data).set_crs("EPSG:4326")

        geometry = [shapely.geometry.Point(xy) for xy in zip(nodes.lon, nodes.lat)]
        nodes = gpd.GeoDataFrame(nodes, crs="EPSG:4326", geometry=geometry)
        nuts_to_nodes = gpd.sjoin(nodes, nuts_data, how='left', op='within')       
        
        # Non Grid Countries get the average availability of all NUTS3 Areas
        availability = pd.DataFrame(index=wind_availability.index)

        non_grid_zones = [z for z in self.zones.index if z not in self.settings["grid_zones"]]
        for zone in non_grid_zones:
            condition_zone = plants.zone == zone
            nuts_areas = nuts_data.loc[nuts_data.country == zone, "name_short"].values
            avg_availability_pv = pv_availability[nuts_areas].mean(axis=1)
            for plant in plants[(plants.technology == "solar") & condition_zone].index:
                availability[plant] = avg_availability_pv
            
            avg_availability_wind = wind_availability[nuts_areas].mean(axis=1)
            for plant in plants[(plants.technology == "wind onshore") & condition_zone].index:
                availability[plant] = avg_availability_wind
        
        for zone in self.settings["grid_zones"]:
            condition_zone = plants.zone == zone
            for plant in plants[(plants.technology == "solar") & condition_zone].index:
                availability[plant] = pv_availability[nuts_to_nodes.loc[plants.loc[plant, "node"], "name_short"]]
            for plant in plants[(plants.technology == "wind onshore") & condition_zone].index:
                availability[plant] = wind_availability[nuts_to_nodes.loc[plants.loc[plant, "node"], "name_short"]]
        
        self.availability = add_timesteps(availability)

    def process_offshore_plants(self):
        
        offshore_availability, offshore_plants = existing_offshore_wind_capacities(self.wdir, self.nodes)
        
        offshore_availability = offshore_availability[(offshore_availability.index >= self.time_horizon["start"]) & \
                                                      (offshore_availability.index < self.time_horizon["end"])]
        
        offshore_availability = add_timesteps(offshore_availability)
        offshore_plants = offshore_plants.drop(["eez_id", "eez_name", "geometry"], axis=1)
        self.plants = pd.concat([self.plants, offshore_plants])
        self.availability = pd.concat([self.availability, offshore_availability], axis=1)

    def create_basic_ntcs(self):
        # from pyhsical cross border flows
        
        pcbf = pd.read_csv(self.wdir.joinpath('data_out/exchange/physical_crossborder_flow.csv'), index_col=0)
        pcbf.utc_timestep = pd.to_datetime(pcbf.utc_timestep).astype('datetime64[ns]')

        self.ntc = pd.DataFrame(index=pd.MultiIndex.from_tuples([(f,t) for (f,t) in itertools.permutations(list(self.zones.index), 2)]))
        max_pcbf = pcbf.groupby(["from_zone", "to_zone"]).max().reset_index()
        self.ntc["ntc"] = 0
        for (f,t) in self.ntc.index:
            self.ntc.loc[(f,t), "ntc"] = max_pcbf.loc[(max_pcbf.from_zone == f) & \
                                                      (max_pcbf.to_zone == t), "value"].max()
        self.ntc = self.ntc.reset_index().fillna(0)
        self.ntc.columns = ["zone_i", "zone_j", "ntc"]
        # Set NTC to zweo if no physical connection exists. 
        for (f,t) in zip(self.ntc.zone_i, self.ntc.zone_j):
            lines, dclines = [], []
            lines += list(self.lines.index[(self.lines.node_i.isin(self.nodes.index[self.nodes.zone == f]))&(self.lines.node_j.isin(self.nodes.index[self.nodes.zone == t])) ])
            lines += list(self.lines.index[(self.lines.node_i.isin(self.nodes.index[self.nodes.zone == t]))&(self.lines.node_j.isin(self.nodes.index[self.nodes.zone == f]))])
            dclines += list(self.dclines.index[(self.dclines.node_i.isin(self.nodes.index[self.nodes.zone == f]))&(self.dclines.node_j.isin(self.nodes.index[self.nodes.zone == t])) ])
            dclines += list(self.dclines.index[(self.dclines.node_i.isin(self.nodes.index[self.nodes.zone == t]))&(self.dclines.node_j.isin(self.nodes.index[self.nodes.zone == f]))])
            
            if len(lines) == 0 and len(dclines) == 0:
                self.ntc.loc[(self.ntc.zone_i == f) & (self.ntc.zone_j == t), "ntc"] = 0
            else: #all(self.ntc.loc[(self.ntc.zone_i == f) & (self.ntc.zone_j == t), "ntc"] == 0):
                self.ntc.loc[(self.ntc.zone_i == f) & (self.ntc.zone_j == t), "ntc"] = self.dclines.loc[dclines, "capacity"].sum()
            
            
        # t, f = "DE", "SE" 
        # tmp_ntc = self.ntc.copy()
        
        # tmp_ntc[(tmp_ntc.zone_i == "DE")&(tmp_ntc.zone_j == "NO")]
        # self.ntc[(tmp_ntc.zone_i == "DE")&(self.ntc.zone_j == "NO")]
        
        # self.dclines        
        
    def connect_small_subnetworks(self):
        """Connect small (<10) node subnetworks to the main network."""
        
        # Create incidence matrix
        A = np.zeros((len(self.lines), len(self.nodes)))
        for i, elem in enumerate(self.lines.index):
            A[i, self.nodes.index.get_loc(self.lines.node_i[elem])] = 1
            A[i, self.nodes.index.get_loc(self.lines.node_j[elem])] = -1
        
        # Find Network Components, i.e. number of unconnectd subnetorks. 
        network_componends = sparse.csgraph.connected_components(np.dot(A.T, A))
        nodes_per_component = {i : len(self.nodes.index[network_componends[1] == i]) for i in range(0, network_componends[0])}
        nodes_in_main_networks = np.isin(network_componends[1], [i for i in nodes_per_component if nodes_per_component[i] > len(self.nodes.index)/10])
        
        # Create synthetic lines between small subnetork and the main network
        add_lines = []
        for i in [subnetwork for subnetwork in nodes_per_component if (1 < nodes_per_component[subnetwork] < 10)]:
            condition_subnetwork = network_componends[1] == i
            tmp_nodes = self.nodes[condition_subnetwork].copy()
            tmp_nodes[["closest_node", "distance_closest_node"]] = "", 0
        
            for n in tmp_nodes.index:
                d = distance(self.nodes.loc[nodes_in_main_networks, "lat"].values.copy(), 
                             self.nodes.loc[nodes_in_main_networks, "lon"].values.copy(), 
                             tmp_nodes.loc[n, "lat"],  tmp_nodes.loc[n, "lon"])
            
                tmp_nodes.loc[n, "closest_node"] = self.nodes[nodes_in_main_networks].index[np.argmin(d)]
                tmp_nodes.loc[n, "distance_closest_node"] = min(d)
        
            node_i = tmp_nodes.index[tmp_nodes.distance_closest_node.argmin()]
            x = (0.301 * 1e-3 * tmp_nodes.loc[node_i, "distance_closest_node"])
            r = (0.06 * 1e-3 * tmp_nodes.loc[node_i, "distance_closest_node"])
            x_pu, r_pu = x / 220**2, r / 220**2
            add_lines.append([node_i, tmp_nodes.loc[node_i, "closest_node"], 
                              tmp_nodes.loc[node_i, "name"], 
                              self.nodes.loc[tmp_nodes.loc[node_i, "closest_node"], "name"],
                              220, tmp_nodes.loc[node_i, "distance_closest_node"], 
                              False, np.nan, "ac_synth", x, r, x_pu, r_pu, 491.556, "online", 1900, False, 1])

        tmp_df = pd.DataFrame(index=["sl" + str(i) for i in range(0, len(add_lines))], 
                                     columns=self.lines.columns, data=add_lines)
        self.lines = self.lines.append(tmp_df)
    
    def data_structure_sheet(self):
        tmp = []
        for d in self.data_structure:
            for (k,v) in self.data_structure[d]["attributes"].items():
                tmp.append([d, k, v, True])
            for (k,v) in self.data_structure[d]["optional attributes"].items():
                tmp.append([d, k, v, False])
        
        return pd.DataFrame(tmp, columns=["data", "attributes", "type", "optional"]).set_index("data")

    def save_to_csv(self, foldername, path=None):
        
        data_structure = self.data_structure_sheet()
        if not path:
            path = self.wdir.joinpath("pomato_datasets").joinpath(foldername)
        else:
            path = path.joinpath(foldername)
        if not path.is_dir():
            path.mkdir()

        print("Saving Data as csv to zip file")
        data_structure.to_csv(path.joinpath("data_structure.csv"))
        for elm in data_structure.index.unique():
            if len(getattr(self, elm)) > 5e4:
                print("Flattening %s", elm)
                cols = getattr(self, elm).columns
                getattr(self, elm).pivot(index=cols[0], columns=cols[1], values=cols[2]).to_csv(path.joinpath(elm + ".csv"))
            else:
                getattr(self, elm).to_csv(path.joinpath(elm + ".csv"))
    
        print("zipping....")
        shutil.make_archive(path, 'zip', path)
        shutil.rmtree(path, ignore_errors=True)
        print("saved!")
        
    def add_dcline(self, node_i, node_j, capacity):
        geod = pyproj.Geod(ellps="WGS84")
        geometry = LineString([Point(self.nodes.loc[node_i, ["lon", "lat"]]), 
                               Point(self.nodes.loc[node_j, ["lon", "lat"]])])
        length = geod.geometry_length(geometry)
        add_dclines = pd.DataFrame(index=["dc" + node_i + "_" + node_j], 
                                   columns=self.dclines.columns,
                                   data=[[node_i, node_j, node_i, node_j, 500, 
                                         length, False, geometry.wkt, "dc", 
                                         capacity, "online", 1900]])
        self.dclines = pd.concat([self.dclines, add_dclines])
        
        
# %%

settings = {
    "grid_zones": ["DE", "FR", "BE", "LU", "NL"],
    # "grid_zones": ["DE"],
    "year": 2020,
    "co2_price": 30,
    "time_horizon": "01.01.2020 - 31.1.2020",
    }

wdir = Path(r"C:\Users\riw\Documents\repositories\pomato_data")
data = PomatoData(wdir, settings)

data.add_dcline("nNO", "nSE", 2000)
data.create_basic_ntcs()
data.plants = data.plants[data.plants.g_max > 5]
data.plants.loc[data.plants.plant_type.isin(["hydro_res", "hydro_psp"]),
        "storage_capacity"] = data.plants.g_max * 24*2

# # data.plants
foldername = "CWE_2030"
# # foldername = "DE_2030"
data.save_to_csv(foldername)

# availability = data.availability
# demand_el = data.demand_el
# dclines = data.dclines
# lines = data.lines
# nodes = data.nodes
# plants = data.plants
# zones = data.zones
# ntc = data.ntc
# technology = data.technology

# plants[plants.node.isna()]

