
import pandas as pd
import numpy as np
import geopandas as gpd




def get_countries_regions_ffe():
    # Download the region types
    url = "http://opendata.ffe.de:3000/rpc/map_region_type?idregiontype=35"
    countries = gpd.read_file(url)
    countries = countries[['name', 'name_short', 'area_m2', "geometry"]].set_index("name_short")
    
    url = "http://opendata.ffe.de:3000/rpc/map_region_type?idregiontype=38"
    regions = gpd.read_file(url)
    regions = regions[['id_region', 'name', 'name_short', 'area_m2', "geometry"]]
    regions["country"] = regions.name_short.str[:2]
    return countries, regions
