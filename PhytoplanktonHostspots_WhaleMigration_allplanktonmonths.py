# %% [markdown]
# ## Wildlife Tracker for Oceans: Real-time map assessment of marine fauna habitat with Phytoplankton hostpots
# 
# This notebook aims to share the algorithm which is used to display the *Phytplankton hotspots* over the marine fauna migration data. Note that "Wildlife Tracker for Oceans" is a geo-framework that must be constrained to any marine fauna tracking data. Thus, the Phytoplankton hotpots are adaptative to the location of any marine fauna tracking data that users may include.
# 
# **Inputs**
# - Great Whales summer migration dataset - Sample [1] `data/azores_whales_records.gpkg`
# - Phytoplankton concentration BlueCloud dataset [2] `phytoplankton/...`
# 
# **Outputs**
# - Map of adaptative hotpots overlapping the marine fauna migration data `output/..`
#     
# ### Some relevant links:
# - Online demo of ["Wildlife Tracker for Oceans"](https://share.streamlit.io/gis4-wildlife/wildlife-tracker-oceans-v0.2pro/main/gis4-oceans.py) If geo-framework is broken contact me to fix it and make it avaible for use: bryanvallejo16@gmail.com
# - 12 month Gallery of [Phytoplankton hotspots](https://gis4-wildlife.github.io/PhytoplanktonGallery-BlueCloud/)
# - Video tutorial about [how to use the geo-framework](https://www.youtube.com/watch?v=IYN5dCJg6os)
# 
# ### Run this notebook online with Binder
# [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gis4-wildlife/PhytoplanktonNotebook-BlueCloud/HEAD)
# 
# ### Some views from "Wildlife Tracker for Oceans"
# 
# #### 1) Hotspots of Phytoplankton concentration
# January
# 
# ![jan](gif/january_hotspot.gif)
# 
# #### 2) Yearly Hotspots of Phytoplankton concentration
# Yearly animation and monthly level
# 
# ![jan](gif/yearly_phytoplankton.gif)
# 
# **References**
# - [1] Silva et al (2014). [Data access](https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study72289508)
# - [2] Sauzede el at (2015). [Blue Cloud Data Access](https://www.blue-cloud.org/demonstrators/zoo-and-phytoplankton-eov-products)
# - [3] WDPA (2022) [Data Access](https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA)

# %% [markdown]
# ______

# %%
# %%
import xarray as xr
import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPoint, Polygon
from keplergl import KeplerGl

# Define depth from 36 levels: [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000]
depth = 0

# Concentration of hotspot from range 0 to 1. Example 0.8 takes over 80% quantile.
hotspot_concentration = 0.8

# %% [markdown]
# **FUNCTIONS**

# %%
def get_WildlifeData():
    '''Functions that gives back a data sample of Great Whale migration
    Return: <geodataframe> crs wgs84    
    '''
    filepath = r'input/azores_whales_records.gpkg'
    geodata = gpd.read_file(filepath, driver='GPKG')
    return geodata

# %%
def get_PhytoHotspot(geodata, depth, hotspot_concentration, zarr_filepath='phytoplankton/combined.zarr'):
    '''Function that gives back CHL concentration for all months, at specific depth, and hotspot concentration
    Input - geodata <geodataframe>, depth: <int>, hotspot_concentration <float>, zarr_filepath <str>
    Output - chl_data: <dataframe>
    '''
    # 1) Bounding box of animal tracking data
    bbox = MultiPoint([point for point in geodata['geometry']]).bounds
    bbox_geom = Polygon([(bbox[0], bbox[1]), (bbox[2], bbox[1]), (bbox[2], bbox[3]), (bbox[0], bbox[3])]).buffer(0.5)
    bbox_bounds = bbox_geom.bounds

    # 2) Reading phytoplankton dataset and filtering
    DS = xr.open_zarr(zarr_filepath).sel(depth=depth)
    DS = DS.sel(longitude=slice(bbox_bounds[0]-10, bbox_bounds[2]+10)).sel(latitude=slice(bbox_bounds[1], bbox_bounds[3]))['soca_chla']
    
    chl_df_list = []
    for month in range(1, 13):
        month_str = ('0000' + str(month))[-2:]
        month_data = DS.sel(time=f'2018-{month_str}').to_dataframe().reset_index(drop=False).dropna()
        month_data['time'] = pd.to_datetime(month_data['time'])  
        percentile = month_data['soca_chla'].quantile(hotspot_concentration)
        month_data = month_data.loc[month_data['soca_chla'] >= percentile]
        month_data = month_data.rename(columns={'soca_chla': 'CHLa_mg_L'})
        chl_df_list.append(month_data)
    
    chl_df = pd.concat(chl_df_list, ignore_index=True)
    return chl_df

# %%
def plot_PhytoMap(geodata, phytoplankton):
    '''Function that returns a map of phytoplankton hotspots and marine wildlife data
    return <mapgl instance>
    '''
    # Ensure time data is in string format
    geodata['timestamp'] = geodata['timestamp'].astype(str)
    phytoplankton['time'] = pd.to_datetime(phytoplankton['time'])
    
    # Extract distinct months
    phytoplankton['month'] = phytoplankton['time'].dt.to_period('M')
    time_min = phytoplankton['time'].min()
    time_max = phytoplankton['time'].max()
    phytoplankton['time'] = phytoplankton['time'].astype(str)
    phytoplankton['month'] = phytoplankton['month'].astype(str)
    # Convert min and max time to strings
    time_min_str = time_min.strftime('%Y-%m-%dT%H:%M:%SZ')
    time_max_str = time_max.strftime('%Y-%m-%dT%H:%M:%SZ')

    config = {
        "version": "v1",
        "config": {
            "visState": {
                "filters": [
                    {
                        "dataId": ["Phytoplankton"],
                        "id": "time_filter",
                        "name": ["time"],
                        "type": "timeRange",
                        "value": [time_min_str, time_max_str],
                        "enlarged": True,
                        "plotType": "histogram",
                        "animationWindow": "free",
                        "yAxis": None
                    }
                ],
                "layers": [
                    {
                        "id": "whales_layer",
                        "type": "geojson",
                        "config": {
                            "dataId": "Great Whales",
                            "label": "Great Whales",
                            "color": [18, 147, 154],
                            "columns": {"geojson": "geometry"},
                            "isVisible": True,
                            "visConfig": {
                                "opacity": 0.8,
                                "strokeOpacity": 0.8,
                                "thickness": 0.5,
                                "colorRange": {
                                    "name": "ColorBrewer Paired-8",
                                    "type": "qualitative",
                                    "category": "ColorBrewer",
                                    "colors": [
                                        "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                                        "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00"
                                    ]
                                },
                                "radius": 10
                            }
                        },
                        "visualChannels": {
                            "colorField": {"name": "wild_id", "type": "string"},
                            "colorScale": "ordinal"
                        }
                    },
                    {
                        "id": "phyto_layer",
                        "type": "grid",
                        "config": {
                            "dataId": "Phytoplankton",
                            "label": "Phytoplankton",
                            "color": [221, 178, 124],
                            "columns": {"lat": "latitude", "lng": "longitude"},
                            "isVisible": True,
                            "visConfig": {
                                "opacity": 0.6,
                                "worldUnitSize": 30,
                                "colorRange": {
                                    "name": "ColorBrewer Greens-9",
                                    "type": "singlehue",
                                    "category": "ColorBrewer",
                                    "colors": [
                                        "#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b",
                                        "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b"
                                    ]
                                },
                                "coverage": 1,
                                "sizeRange": [0, 500],
                                "elevationScale": 5,
                                "colorAggregation": "average",
                                "sizeAggregation": "count",
                                "enable3d": False
                            }
                        },
                        "visualChannels": {
                            "colorField": {"name": "CHLa_mg_L", "type": "real"},
                            "colorScale": "quantile"
                        }
                    }
                ],
                "interactionConfig": {
                    "tooltip": {
                        "fieldsToShow": {
                            "Great Whales": [
                                {"name": "timestamp", "format": None},
                                {"name": "location_long", "format": None},
                                {"name": "individual_id", "format": None},
                                {"name": "tag_id", "format": None},
                                {"name": "wild_id", "format": None}
                            ],
                            "Phytoplankton": [
                                {"name": "time", "format": None},
                                {"name": "depth", "format": None},
                                {"name": "CHLa_mg_L", "format": None}
                            ]
                        },
                        "compareMode": False,
                        "enabled": True
                    },
                    "brush": {"enabled": False},
                    "geocoder": {"enabled": False},
                    "coordinate": {"enabled": False}
                },
                "layerBlending": "normal",
                "splitMaps": [],
                "animationConfig": {"currentTime": None, "speed": 1}
            },
            "mapState": {
                "bearing": 0,
                "dragRotate": False,
                "latitude": 39.08282277285858,
                "longitude": -47.56469878817516,
                "pitch": 0,
                "zoom": 1.7759365278629409,
                "isSplit": False
            },
            "mapStyle": {
                "styleType": "satellite",
                "visibleLayerGroups": {
                    "label": True,
                    "road": True,
                    "building": True,
                    "water": True,
                    "land": True
                },
                "threeDBuildingColor": [9.665468314072013, 17.18305478057247, 31.1442867897876],
                "mapStyles": {}
            }
        }
    }

    # Create the map
    Map = KeplerGl(height=800)
    Map.add_data(geodata, 'Great Whales')
    Map.add_data(phytoplankton, 'Phytoplankton')
    Map.config = config

    return Map

# Main execution
geodata = get_WildlifeData()
phyto_data = get_PhytoHotspot(geodata, depth, hotspot_concentration)
map_instance = plot_PhytoMap(geodata, phyto_data)
map_instance.save_to_html(file_name='phyto_map_all_months.html')
