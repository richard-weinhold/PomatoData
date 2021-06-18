"""Preprocess GridKit Data.

Using the pyPSA data from https://github.com/PyPSA/pypsa-eur/tree/master/data/entsoegridkit

"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import time
import geopy
import geopandas as gpd
from shapely.geometry import Point
from geopy.exc import GeocoderQueryError, GeocoderTimedOut, GeocoderQuotaExceeded, GeocoderServiceError

def return_symbol_from_tag(tags, symbol):
    # tags = nodes.tags
    # symbol = "asd"
    df_return = pd.Series(index=tags.index, dtype=float)
    tmp = tags.str.split(",", expand=True)
    for col in tmp.columns:
        df_return[tmp.loc[:, col].str.contains(symbol, na=False)] = tmp.loc[:, col][tmp.loc[:, col].str.contains(symbol, na=False)]
    
    if df_return.isna().all():
        return df_return.values
    else:    
        return df_return.str.split("=>", expand=True).iloc[:, 1].str.replace('"', '').values

def distance(lat_nodes, lon_nodes, lat_plants, lon_plants):
    """Return vector of distances in km from plant to all nodes in zone."""
    R = 6373.0
    lat_nodes *= np.pi / 180
    lat_plants *= np.pi / 180
    lon_plants *= np.pi / 180
    dlon = lon_nodes - lon_plants
    dlat = lat_nodes - lat_plants
    a = np.sin(dlat / 2)**2 + np.cos(lat_nodes) * np.cos(lat_plants) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R * c

def check_all_nodes_in_shape(wdir, nodes, zone):
    country_shapes = gpd.read_file(str(wdir.joinpath('data/demand_node_data/NUTS_2013_01M_SH/data/NUTS_RG_01M_2013.shp')))
    country_shape = country_shapes.loc[country_shapes['STAT_LEVL_'] == 0]
    
    tmp_shape = country_shape.loc[country_shape['NUTS_ID'] == zone, "geometry"].values[0]
    tmp_points = [Point(xy) for xy in zip(nodes.loc[nodes.zone ==  zone, "lon"], nodes.loc[nodes.zone ==  zone, "lat"])]
    tmp_contrains = [not p.within(tmp_shape) for p in tmp_points]
    tmp = nodes[(nodes.zone ==  zone)].loc[tmp_contrains].loc[nodes.demand]
    if not tmp.empty: 
        print("Nodes not in zone", " ".join(tmp.index))
    return nodes[(nodes.zone ==  zone)].loc[tmp_contrains]

# %%
def process_gridkit_data(gridkit_filepath, version="jan_2020"):
    # gridkit_filepath = Path(r"C:\Users\riw\Documents\repositories\pomato_data\data_in\GridKit")
    nodes = pd.read_csv(gridkit_filepath.joinpath(version + "/buses.csv"),
                        quotechar="'", true_values='t', false_values='f')
    
    lines = pd.read_csv(gridkit_filepath.joinpath(version + "/lines.csv"),
                        quotechar="'", true_values='t', false_values='f')
    
    transformers = pd.read_csv(gridkit_filepath.joinpath(version + "/transformers.csv"),
                               quotechar="'", true_values='t', false_values='f', index_col=0)
    
    dclines = pd.read_csv(gridkit_filepath.joinpath(version + "/links.csv"),
                          quotechar="'", true_values='t', false_values='f')
    
    nodes["zone"] = return_symbol_from_tag(nodes.tags, "country")
    nodes["name"] = return_symbol_from_tag(nodes.tags, "name_eng")
    nodes = nodes.set_index("bus_id")
    
    lines["from_zone"] = return_symbol_from_tag(lines.tags, "country_1")
    lines["to_zone"] = return_symbol_from_tag(lines.tags, "country_2")
    
    for node in nodes.index[nodes.voltage.isna()]:
        tmp_lines = lines[(lines.bus0 == node) | (lines.bus1 == node)][["bus0", "bus1"]].copy()
        if len(tmp_lines) == 0:
            nodes.loc[node, "voltage"] = 500
        elif len(tmp_lines) ==  1:
            nodes.loc[node, "voltage"] = tmp_lines.voltage.values[0]
        else:
            print(f"Node {node} is inbetween {' '.join(list(tmp_lines.voltages.values))}. Choosing first")
            nodes.loc[node, "voltage"] = tmp_lines.voltage.values[0]
    
    
    dclines["name"] = return_symbol_from_tag(dclines.tags, "text_")
    dclines["from_zone"] = return_symbol_from_tag(dclines.tags, "country_1")
    dclines["to_zone"] = return_symbol_from_tag(dclines.tags, "country_2")
    
    lines["name_i"] = nodes.name[lines.bus0].values
    lines["name_j"] = nodes.name[lines.bus1].values
    transformers["name_i"] = nodes.name[transformers.bus0].values
    transformers["name_j"] = nodes.name[transformers.bus1].values
    
    dclines["name_i"] = nodes.name[dclines.bus0].values
    dclines["name_j"] = nodes.name[dclines.bus1].values
    
    transformers["voltage"] = ""
    transformers["voltage"] = nodes.voltage[transformers.bus0.values].astype(int).map(str).values + "/" + nodes.voltage[transformers.bus1.values].astype(int).map(str).values
    
    cables = lines[lines.underground].copy()
    
    # %%
    
    nodes = nodes[["station_id", "voltage", "name", "y", "x", "zone", "symbol"]]
    nodes.columns = ['substation', 'voltage', 'name', 'lat', 'lon', 'zone', 'info', ]
    
    lines = lines[["line_id", "bus0", "bus1", "name_i", "name_j", 'voltage',
                   'circuits', 'length', "under_construction", "geometry"]].set_index("line_id")
    
    dclines = dclines[["link_id", "bus0", "bus1", "name_i", "name_j",
                       'length', "under_construction", "geometry"]].set_index("link_id")
    
    transformers.columns = ["node_i", "node_j", "name_i", "name_j", 'voltage']
    
    lines.columns = ["node_i", "node_j", "name_i", "name_j",
                     'voltage', 'circuits', 'length', "under_construction", "geometry"]
    
    dclines.columns = ["node_i", "node_j", "name_i", "name_j", 'length', "under_construction", "geometry"]
    
    dclines["circuits"] = ""
    dclines["voltage"] = ""
    dclines["x"] = ""
    dclines["r"] = ""
    dclines["x_pu"] = ""
    dclines["r_pu"] = ""
    dclines["capacity"] = 1000
    dclines["technology"] = "dc"
    dclines["status"] = dclines["under_construction"].map({True: "planned", False: "online"})
    dclines["commissioned"] = dclines["status"].map({"online": 1900, "planned": 2050})
    
    # https://github.com/PyPSA/EnergyModels.jl/blob/master/src/data/pypsa.lines.types.csv 
    line_type = {132: {"r_per_km": 0.0949, "x_per_km": 0.38, "i_nom": 0.74},  # 305-AL1/39-ST1A 110.0
                 220: {"r_per_km": 0.06, "x_per_km": 0.301, "i_nom": 1.29},   # Al/St 240/40 2-bundle 220.0
                 300: {"r_per_km": 0.04, "x_per_km": 0.265, "i_nom": 1.935},  # Al/St 240/40 3-bundle 300.0
                 380: {"r_per_km": 0.03, "x_per_km": 0.246, "i_nom": 2.58},   # Al/St 240/40 4-bundle 380.0
                 500: {"r_per_km": 0.02, "x_per_km": 0.222, "i_nom": 3.225},  # Hypothetical 5-bundle 300.0 * 5/3
                 750: {"r_per_km": 0.01, "x_per_km": 0.202, "i_nom": 3.87}}   # Hypothetical 6-bundle 300.0 * 6/3
    
    cable_type = {132: {"r_per_km": 0.0477, "x_per_km": 0.1298, "i_nom": 0.588},    # Oeding and Oswald
                  220: {"r_per_km": 0.0344, "x_per_km": 0.1261, "i_nom": 1.29},     # Oeding and Oswald
                  380: {"r_per_km": 0.0312, "x_per_km": 0.1281, "i_nom": 1.935}}    # Oeding and Oswald
    
    lines["technology"] = "ac_line"
    lines.loc[cables.line_id, "technology"] = "ac_cable"
    lines["x"] = lines.voltage.map({k: line_type[k]["x_per_km"] for k in line_type}) * lines.length * 1e-3
    lines["r"] = lines.voltage.map({k: line_type[k]["r_per_km"] for k in line_type}) * lines.length * 1e-3
    lines["x_pu"] = lines["x"]/lines["voltage"]**2
    lines["r_pu"] = lines["r"]/lines["voltage"]**2
    
    # lines.loc[cables.line_id, "x"] = lines.loc[cables.line_id].voltage.map({k: line_type[k]["x_per_km"] for k in cable_type})
    # lines.loc[cables.line_id, "r"] = lines.loc[cables.line_id].voltage.map({k: line_type[k]["r_per_km"] for k in cable_type})
    
    lines["capacity"] = np.sqrt(3) * lines.voltage * lines.voltage.map({k: line_type[k]["i_nom"] for k in line_type})
    lines["status"] = lines["under_construction"].map({True: "planned", False: "online"})
    lines["commissioned"] = lines["status"].map({"online": 1900, "planned": 2050})
    
    
    # https://circuitglobe.com/what-is-resistance-and-reactance-of-the-transformer.html
    # PyPSA Code
    # t["r"] = t["vscr"] /100.
    # t["x"] = np.sqrt((t["vsc"]/100.)**2 - t["r"]**2)
    
    transformers.voltage.unique()
    # https://pandapower.readthedocs.io/en/v1.3.0/std_types/basic.html#transformers
    # https://github.com/lthurner/pandapower/blob/v1.2.2/pandapower/std_types.py
    transformer_type = {
        "110/220": {"Un": 220, "vsc": 12, "vscr": 0.26, "Sn": 100},   # 100 MVA 220/110 kV
        "132/220": {"Un": 220, "vsc": 12, "vscr": 0.26, "Sn": 100},   # 100 MVA 220/110 kV
        "132/300": {"Un": 300, "vsc": 12.2, "vscr": 0.25, "Sn": 160}, # 160 MVA 380/110 kV
        "220/300": {"Un": 300, "vsc": 18.5, "vscr": 0.25, "Sn": 300}, # simbench Typ_x_380/220
        "110/380": {"Un": 380, "vsc": 12.2, "vscr": 0.25, "Sn": 300}, # 160 MVA 380/110 kV
        "132/380": {"Un": 380, "vsc": 12.2, "vscr": 0.25, "Sn": 300}, # 160 MVA 380/110 kV
        "220/380": {"Un": 380, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        "300/380": {"Un": 380, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
    
        "132/500": {"Un": 500, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        "220/500": {"Un": 500, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        "300/500": {"Un": 500, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        "380/500": {"Un": 500, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        "132/750": {"Un": 750, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        "220/750": {"Un": 750, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        "300/750": {"Un": 750, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        "380/750": {"Un": 750, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        "500/750": {"Un": 750, "vsc": 18.5, "vscr": 0.25, "Sn": 600}, # simbench Typ_x_380/220
        }
    
    # Un = 380e3
    # Sn = 400e6
    # ur = 0.25/100
    # uk = 18.5/100
    
    # x = (uk**2 - ur**2)**0.5 / Sn * Un**2
    # r = ur / Sn * Un**2
    
    # x = uk * (Un/Ub)**2 * 1/Sn
    
    transformers["circuits"] = ""
    transformers["length"] = 1
    transformers["geometry"] = ""
    
    transformers["ur"] = transformers.voltage.map({k: transformer_type[k]["vscr"] for k in transformer_type}) / 100
    transformers["uk"] = transformers.voltage.map({k: transformer_type[k]["vsc"] for k in transformer_type}) / 100
    transformers["Un"] = transformers.voltage.map({k: transformer_type[k]["Un"] for k in transformer_type}) * 1e3
    transformers["Sn"] = transformers.voltage.map({k: transformer_type[k]["Sn"] for k in transformer_type}) * 1e6
    
    transformers["r"] = (transformers["ur"] * transformers["Un"]**2) / transformers["Sn"]
    transformers["z"] = (transformers["uk"] * transformers["Un"]**2) / transformers["Sn"]
    transformers["x"] = np.sqrt(transformers["z"]**2 - transformers["r"]**2)
    
    transformers["x_pu"] = transformers["x"] / transformers["Sn"]
    transformers["r_pu"] = transformers["r"] / transformers["Sn"]
    
    transformers["capacity"] = transformers["Sn"]
    transformers["technology"] = "transformer"
    transformers["status"] = "online"
    transformers["commissioned"] = 1900
    transformers["under_construction"] = ""
    
    nodes["demand"] = True
    nodes.loc[nodes["info"] != "Substation", "demand"] = False
    
    for node in nodes.index[nodes.zone.isna()]:
        tmp_lines = lines[(lines.node_i == node) | (lines.node_j == node)][["node_i", "node_j"]].copy()
        tmp_zones = list(nodes.loc[tmp_lines.node_i, "zone"]) + list(nodes.loc[tmp_lines.node_j, "zone"])
        zones = [zone for zone in set(tmp_zones) if zone not in [" ", None, np.nan]]
        if len(zones) == 1 or len(set(zones)) == 1:
            nodes.loc[node, "zone"] = zones[0]
        else:
            nodes.loc[node, "zone"] = np.nan
    
    ## Make sure index type is String and not int
    lines.index = lines.index.astype(str)
    lines.index = "l" + lines.index
    dclines.index = dclines.index.astype(str)
    dclines.index = "dc" + dclines.index
    transformers.index = transformers.index.astype(str)
    transformers.index = "tr" + transformers.index
    
    lines.node_i = lines.node_i.astype(str)
    lines.node_j = lines.node_j.astype(str)
    lines.node_i = "n" + lines.node_i
    lines.node_j = "n" + lines.node_j
    
    transformers.node_i = transformers.node_i.astype(str)
    transformers.node_j = transformers.node_j.astype(str)
    transformers.node_i = "n" + transformers.node_i
    transformers.node_j = "n" + transformers.node_j
    
    dclines.node_i = dclines.node_i.astype(str)
    dclines.node_j = dclines.node_j.astype(str)
    dclines.node_i = "n" + dclines.node_i
    dclines.node_j = "n" + dclines.node_j
    
    nodes.index = nodes.index.astype(str)
    nodes.index = "n" + nodes.index
    
    nodes.substation = "s" + nodes.substation.astype(str)
    
    nodes.index.name = 'index'
    lines.index.name = 'index'
    dclines.index.name = 'index'
    transformers.index.name = 'index'
    
    nodes.index[nodes.zone.isna()]
    
    # %%
    print("Number of nodes without zone:", len(nodes.loc[(nodes.zone.isna())|(nodes.zone == " ")]))
    
    from geopy.geocoders import Nominatim#, ArcGIS
    geolocator = Nominatim(user_agent="TUBERLIN WIP")
    # location = geolocator.reverse((55.328611, 12.292824))
    
    for counter, node in enumerate(nodes.loc[(nodes.zone.isna())|(nodes.zone == " ")].index):
        if counter%10 == 0:
            print(counter)
        while True:
            try:
                location = geolocator.reverse((nodes.loc[node, "lat"], nodes.loc[node, "lon"]))
                nodes.loc[node, "zone"] = str.upper(location.raw["address"]["country_code"])
                break
            except GeocoderTimedOut:
                print("Code 5: Time Out - Sleep")
                time.sleep(2)
            except GeocoderServiceError:
                print("Code 6: Service Error - Sleep")
                time.sleep(2)
            except KeyboardInterrupt:
                print(f"Last address: {location.address}")
                raise
            except:
                print("Unexpected error:", sys.exc_info()[0])
                latlon = (nodes.loc[node, "lat"], nodes.loc[node, "lon"])
                print(f"Coordinates causing problems: {latlon}")
                print(f"Address causing problems: {location.address}")
                print("Zone set to 'None', fix manually")
                nodes.loc[node, "zone"] = "None"
                break
                
    nodes[nodes.zone == ""]
    if version == "jun_2020":
        nodes.loc["n7278", "zone"] = "DK"
    
    print("Number of nodes without zone:", len(nodes.index[nodes.zone.isna()]))
    
    # %%
    
    # nodes.loc["n6771", :]                 
    # check_all_nodes_in_shape(nodes, "LU")
    if version == "jan_2020":
        # Nodes fix DE
        nodes.loc["n3848", ["lat", "lon"]] = 47.593841, 7.874779
        nodes.loc["n5542", ["lat", "lon"]] = 54.143660, 12.128323
        nodes.loc["n5683", ["lat", "lon"]] = 54.334829, 10.118415
        nodes.loc["n5729", ["lat", "lon"]] = 53.6471, 8.0
        
        # Fix FR
        nodes.loc["n2482", ["lat", "lon"]] = 43.368001, 3.539757 
        
        #Fix NL
        nodes.loc["n6771", ["lat", "lon"]] = 53.332232, 6.925237
        
        # Nodes fix NO
        nodes.loc["n6489", ["lat", "lon"]] =  60.458071, 5.274477 
        nodes.loc["n6478", ["lat", "lon"]] =  61.491851, 5.749391 
        
        lines.loc["l548", "circuits"] = 2
        lines.loc["l965", "capacity"] = 500
        
        # Double Line Redwitz Altenfeld
        redwitz_altenfeld = lines.loc[["l769"]].copy()
        redwitz_altenfeld.loc[:, "node_j"] = "n4554"
        redwitz_altenfeld["index"] = ["red_alt"]
        redwitz_altenfeld.set_index("index", drop=True, inplace=True)
        lines = lines.append(redwitz_altenfeld)
        
        # New Circuit Vieselbach - Lauchstädt
        lines.loc["l875", "circuits"] = 3
        # Two Circuits Emsland
        lines.loc["l957", "circuits"] = 2
        
        # Two circuits nuclear PP Gravelines
        lines.loc["l5528", "circuits"] = 3
        
        # 3 Circuits for SE-DE Interconnector 
        lines.loc["l725", "circuits"] = 3
        lines.loc["l897", "circuits"] = 3
        
        # FR-ES Interconnector 
        lines.loc["l5533", "circuits"] = 2
        lines.loc["l4873", "circuits"] = 2
        lines.loc["l4870", "circuits"] = 2

        # Line in France (Nice), wich has bad demand matched. 
        lines.loc["l5373", "circuits"] = 2
        # other line in france
        lines.loc["l5726", "circuits"] = 2
        
        # Paris Ring
        lines.loc["l5421", "circuits"] = 4
        lines.loc["l5421", "capacity"] *= 2
        lines.loc["l5524", "circuits"] = 4
        lines.loc["l5524", "capacity"] *= 2

        lines.loc["l5435", "circuits"] = 4
        lines.loc["l5435", "capacity"] *= 2
        
    if version == "jun_2020":
        nodes.loc["n3513", ["lat", "lon"]] =  47.579165, 7.810306 
        nodes.loc["n5666", ["lat", "lon"]] = 53.639167, 8.035150
        nodes.loc["n5665", ["lat", "lon"]] = 53.639167, 8.035150
        nodes.loc["n5730", ["lat", "lon"]] = 54.360133, 10.142023
    
    # for zone in [zone for zone in nodes.zone.unique() if zone in possible_zones][0:3]:
    # Elba -> IT
    nodes.loc[(nodes.lat < 43)&(nodes.lon > 7.6)&(nodes.zone == "FR"), "zone"] = "IT"
    nodes.loc[:, "zone"] = nodes.zone.replace({"GB": "UK"})
    
    # %%
    def n_points_on_curcle(n, r, center):
        alpha = 2*np.pi/n
        return [(center[0] + r*np.cos(k*alpha), center[1] + r*np.sin(k*alpha)) for k in range(0, n)]
    
    tmp_nodes = nodes.copy()
    tmp_nodes.loc[:, ["lat", "lon"]] = nodes.loc[:,  ["lat", "lon"]].apply(round, ndigits=3)
    tmp_nodes["counter"] = 1
    tmp_nodes = tmp_nodes[["counter", "lat", "lon"]].groupby(["lat", "lon"]).sum().reset_index()
    
    relevant_coords = tmp_nodes.loc[tmp_nodes.counter > 1,  ["lat", "lon"]].apply(tuple, axis=1)
    nodes.loc[:,  ["lat", "lon"]].apply(round, ndigits=3).apply(tuple, axis=1)
    
    i = 0
    for coord in relevant_coords:
        i += 1
        if i%100 == 0:
            print(i)
        cond = nodes[["lat", "lon"]].apply(round, ndigits=3).apply(tuple, axis=1) == coord
        if sum(cond) > 1:
            nodes.loc[cond, ["lat", "lon"]] = n_points_on_curcle(sum(cond), 0.01, coord)
        else:
            print(coord)

    lines_combined = pd.concat([lines, dclines[lines.columns], transformers[lines.columns]], axis=0)
    return nodes, lines_combined     
    

# %%

if __name__ == "__main__":
    
    gridkit_filepath = Path(r"C:/Users/riw/Documents/repositories/pomato_data/data_in/GridKit")
    nodes, lines = process_gridkit_data(gridkit_filepath)
    data_out_folder = Path(r"C:/Users/riw/Documents/repositories/pomato_data/data_out")
    data_in_folder = Path(r"C:/Users/riw/Documents/repositories/pomato_data/data_in")
        
    add_dclines = pd.read_csv(data_in_folder.joinpath("grid/add_dclines.csv"), index_col=0)
    tmp_lines = pd.concat([lines, add_dclines], axis=0)
    tmp_lines.to_csv(data_out_folder.joinpath("lines/lines.csv"))
    nodes.to_csv(data_out_folder.joinpath("nodes/nodes.csv"))
    
    # lines = pd.read_csv(data_out_folder.joinpath("lines/lines.csv"), index_col=0)
    # add_dclines = pd.read_csv(data_in_folder.joinpath("grid/add_dclines.csv"), index_col=0)
    # lines = pd.concat([lines, add_dclines], axis=1)
    # t = check_all_nodes_in_shape(nodes, "NO")
    # t = check_all_nodes_in_shape(nodes, "DE")

    # %% Checking Staff
    
    # cwe = ["FR", "BE", "NL", "LU", "DE"]
    # cwe_nodes = nodes[nodes.zone.isin(cwe)].index
    # de_nodes = nodes[nodes.zone == "DE"].index
    
    # lines["sum"] = 1
    # lines.loc[lines.node_i.isin(de_nodes) | lines.node_j.isin(de_nodes), ["sum", "voltage"]].groupby("voltage").sum()
    # lines.loc[lines.node_i.isin(cwe) | lines.node_j.isin(cwe), ["sum", "voltage"]].groupby("voltage").sum()
    
    # lines[(lines.node_i.isin(cwe_nodes) | lines.node_j.isin(cwe_nodes))]
    # lines[lines.node_i.isin(de_nodes) | lines.node_j.isin(de_nodes) ]
    
    # transformers[transformers.node_i.isin(cwe_nodes) | transformers.node_j.isin(cwe_nodes) ].voltage.unique()
    # transformers[transformers.node_i.isin(de_nodes) | transformers.node_j.isin(de_nodes) ].voltage.unique()
    
    
    # dclines[dclines.node_i.isin(cwe_nodes) | dclines.node_j.isin(cwe_nodes) ]
    # dclines[dclines.node_i.isin(de_nodes) | dclines.node_j.isin(de_nodes) ]
    
    
    # nodes[nodes.name == "Neuenhagen"]
    # lines[lines.node_i.isin(r) | lines.node_j.isin(r)]
    
    # r = ["n4499"]
    # r = ["n4545"]
    # nodes.loc["n7610"]
    # lines.loc["l9602"]
    
    # transformers.loc["tr12371"]
    # transformers.loc["tr12366"]
    
    # transformers[transformers.node_i.isin(r) | transformers.node_j.isin(r) ]
    
    # %%
    
    # for dc in dclines.index:
    #     line_nodes =  list(lines.node_j) + list(lines.node_i)
    #     connected_dclines = dclines[(dclines.node_i == dclines.loc[dc, "node_i"]) | (dclines.node_j == dclines.loc[dc, "node_i"])]
    #     if not (dclines.loc[dc, "node_i"] in line_nodes):
    #         if len(connected_dclines) > 1:
    #             print(f"{dc} not connected to a ac line at node_i, but with dc line")
    #         else:
    #             print(f"{dc} not connected to a ac or dc line at node_i")
    #     if not (dclines.loc[dc, "node_j"] in line_nodes):
    #         if len(connected_dclines) > 1:
    #             print(f"{dc} not connected to a ac line at node_j, but with dc line")
    #         else:
    #             print(f"{dc} not connected to a ac or dc line at node_j")
    # print("done")


