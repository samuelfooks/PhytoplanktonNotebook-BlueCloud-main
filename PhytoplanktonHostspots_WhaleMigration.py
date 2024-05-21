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
import geopandas as gpd
import xarray as xr                                                # pip install xarray
import os
from shapely.geometry import MultiPoint, Polygon, Point
from keplergl import KeplerGl                                      # pip install keplergl

# import cmocean                        # pip install cmocean
# import netCDF4                        # pip install netCDF4

# %% [markdown]
# **PARAMETERS**
# 
# If you change parameters, please run the whole notebook

# %%
# define month from range 1 to 12. Example July = 7
month = 10

# define depth from 36 levels: [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000]
depth = 5

# concentration of hotspot from range 0 to 1. Example 0.8 takes over 80% quantile. Test it with 0 and 0.8 :) 
hotspot_concentration = 0.8

# %% [markdown]
# **FUNCTIONS**

# %%
def get_WildlifeData():
    '''Functions that gives back a data sample of Great Whale migration
    Return: <geodataframe> crs wgs84    
    '''
    
    filepath = r'input/azores_whales_records.gpkg'
    
    geodata = gpd.read_file(filepath, driver = 'GPKG')
    
    return geodata

# %%
def get_PhytoHotspot(geodata, month, depth, hotspot_concentration):
    '''Functions that gives back CHL concentration in specific month, at specific concentration, at specific depth
    Input - geodata <geodataframe>, month: <int>, depth: <int>, hotspot_concentration <float> 
    Output - chl_data: <dataframe>
    '''
    
    # 1) Bounding box of animal tracking data
    
    bbox = MultiPoint([point for point in geodata['geometry']]).bounds
    
    bbox_geom = Polygon([(bbox[0], bbox[1]), (bbox[2], bbox[1]),  (bbox[2], bbox[3]), (bbox[0], bbox[3])]).buffer(0.5)
    bbox_bounds = bbox_geom.bounds

    # 2) Filepath of Phytoplankton concentration
    
    month_str = ('0000' + str(month))[-2:]
    chl_filepath = fr'phytoplankton/{month}/SOCA_CHLA_glo_bgc3d_rep_2018_{month_str}_P20210222.nc'
    
    # 3) Reading phytoplankton dataset and filtering
    
    DS = xr.open_dataset(chl_filepath).sel(depth=depth)
    DS = DS.sel(longitude=slice(bbox_bounds[0]-10,bbox_bounds[2]+10)).sel(latitude=slice(bbox_bounds[1],bbox_bounds[3]))['soca_chla']

    chl_df = DS.to_dataframe()
    chl_df = chl_df.reset_index(drop=False).dropna()
    
    percentile = chl_df['soca_chla'].quantile(hotspot_concentration)
    chl_df = chl_df.loc[chl_df['soca_chla']>=percentile]
        
    chl_df = chl_df.rename(columns = {'soca_chla': 'CHLa_mg_L'})
    
    return chl_df    

# %%
def plot_PhytoMap(geodata, phytoplankton):
    '''Function that return a map of phytoplankton hotspots and marine wildlife data
    return <mapgl instance>
    '''
    config = {
  "version": "v1",
  "config": {
    "visState": {
      "filters": [
        {
          "dataId": [
            "Great Whales"
          ],
          "id": "hflv90e0g",
          "name": [
            "timestamp"
          ],
          "type": "timeRange",
          "value": [
            1211810400000,
            1222923509000.0002
          ],
          "enlarged": True,
          "plotType": "histogram",
          "animationWindow": "free",
          "yAxis": None
        }
      ],
      "layers": [
        {
          "id": "cm7q9ak",
          "type": "geojson",
          "config": {
            "dataId": "Great Whales",
            "label": "Great Whales",
            "color": [
              18,
              147,
              154
            ],
            "columns": {
              "geojson": "geometry"
            },
            "isVisible": True,
            "visConfig": {
              "opacity": 0.8,
              "strokeOpacity": 0.8,
              "thickness": 0.5,
              "strokeColor": None,
              "colorRange": {
                "name": "ColorBrewer Paired-8",
                "type": "qualitative",
                "category": "ColorBrewer",
                "colors": [
                  "#a6cee3",
                  "#1f78b4",
                  "#b2df8a",
                  "#33a02c",
                  "#fb9a99",
                  "#e31a1c",
                  "#fdbf6f",
                  "#ff7f00"
                ]
              },
              "strokeColorRange": {
                "name": "Global Warming",
                "type": "sequential",
                "category": "Uber",
                "colors": [
                  "#5A1846",
                  "#900C3F",
                  "#C70039",
                  "#E3611C",
                  "#F1920E",
                  "#FFC300"
                ]
              },
              "radius": 10,
              "sizeRange": [
                0,
                10
              ],
              "radiusRange": [
                0,
                50
              ],
              "heightRange": [
                0,
                500
              ],
              "elevationScale": 5,
              "stroked": False,
              "filled": True,
              "enable3d": False,
              "wireframe": False
            },
            "hidden": False,
            "textLabel": [
              {
                "field": None,
                "color": [
                  255,
                  255,
                  255
                ],
                "size": 18,
                "offset": [
                  0,
                  0
                ],
                "anchor": "start",
                "alignment": "center"
              }
            ]
          },
          "visualChannels": {
            "colorField": {
              "name": "wild_id",
              "type": "string"
            },
            "colorScale": "ordinal",
            "sizeField": None,
            "sizeScale": "linear",
            "strokeColorField": None,
            "strokeColorScale": "quantile",
            "heightField": None,
            "heightScale": "linear",
            "radiusField": None,
            "radiusScale": "linear"
          }
        },
        {
          "id": "lbhpz34",
          "type": "grid",
          "config": {
            "dataId": "Phytoplankton",
            "label": "Point",
            "color": [
              221,
              178,
              124
            ],
            "columns": {
              "lat": "latitude",
              "lng": "longitude"
            },
            "isVisible": True,
            "visConfig": {
              "opacity": 0.6,
              "worldUnitSize": 30,
              "colorRange": {
                "name": "ColorBrewer Greens-9",
                "type": "singlehue",
                "category": "ColorBrewer",
                "colors": [
                  "#f7fcf5",
                  "#e5f5e0",
                  "#c7e9c0",
                  "#a1d99b",
                  "#74c476",
                  "#41ab5d",
                  "#238b45",
                  "#006d2c",
                  "#00441b"
                ]
              },
              "coverage": 1,
              "sizeRange": [
                0,
                500
              ],
              "percentile": [
                0,
                100
              ],
              "elevationPercentile": [
                0,
                100
              ],
              "elevationScale": 5,
              "colorAggregation": "average",
              "sizeAggregation": "count",
              "enable3d": False
            },
            "hidden": False,
            "textLabel": [
              {
                "field": None,
                "color": [
                  255,
                  255,
                  255
                ],
                "size": 18,
                "offset": [
                  0,
                  0
                ],
                "anchor": "start",
                "alignment": "center"
              }
            ]
          },
          "visualChannels": {
            "colorField": {
              "name": "CHLa_mg_L",
              "type": "real"
            },
            "colorScale": "quantile",
            "sizeField": None,
            "sizeScale": "linear"
          }
        }
      ],
      "interactionConfig": {
        "tooltip": {
          "fieldsToShow": {
            "Great Whales": [
              {
                "name": "timestamp",
                "format": None
              },
              {
                "name": "location_long",
                "format": None
              },
              {
                "name": "individual_id",
                "format": None
              },
              {
                "name": "tag_id",
                "format": None
              },
              {
                "name": "wild_id",
                "format": None
              }
            ],
            "Phytoplankton": [
              {
                "name": "time",
                "format": None
              },
              {
                "name": "depth",
                "format": None
              },
              {
                "name": "CHLa_mg_L",
                "format": None
              }
            ]
          },
          "compareMode": False,
          "compareType": "absolute",
          "enabled": True
        },
        "brush": {
          "size": 0.5,
          "enabled": False
        },
        "geocoder": {
          "enabled": False
        },
        "coordinate": {
          "enabled": False
        }
      },
      "layerBlending": "normal",
      "splitMaps": [],
      "animationConfig": {
        "currentTime": None,
        "speed": 1
      }
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
      "topLayerGroups": {},
      "visibleLayerGroups": {
        "label": True,
        "road": True,
        "border": False,
        "building": True,
        "water": True,
        "land": True,
        "3d building": False
      },
      "threeDBuildingColor": [
        9.665468314072013,
        17.18305478057247,
        31.1442867897876
      ],
      "mapStyles": {}
    }
  }
}
    
    
    geodata['timestamp'] = geodata['timestamp'].astype(str)
    phytoplankton['time'] = phytoplankton['time'].astype(str)
    
#     month = str(phytoplankton.time.unique()[0]).split('-')[1]
    
    Map = KeplerGl(height = 800)
    Map.add_data(geodata, 'Great Whales')
    Map.add_data(phytoplankton, 'Phytoplankton')
    
    Map.config = config
    Map.save_to_html(file_name = f'outputhtml/phyto_map_{month}.html')
    return Map

# %%
def save_PhytoMap(month, depth, hotspot_concentration):
    '''Function that saves the current map 
    Return <html> in output folder
    '''
    
    # cwd
    home = os.getcwd()
    output_folder = os.path.join(home, 'output')

    if not os.path.exists(output_folder):
        os.makedirs(output_folder) 

    # try saving
    Map.save_to_html(file_name = f'outputhtml/Phytoplankton_hotspot-{hotspot_concentration}_month-{month}_depth-{depth}_.html')

# %% [markdown]
# #### 1) Reading data
# This demo works with Phytoplankton data products obtained from Blue Cloud Vlabs analyzed for 2018. Products are stored in folder phytoplankton.

# %%
# Marine fauna tracking data sample
geodata = get_WildlifeData()
geodata.head()

# %%
# Phytoplankton hotspots
phytoplankton = get_PhytoHotspot(geodata, month, depth, hotspot_concentration)
phytoplankton.head()

# %% [markdown]
# #### 2) Ploting data 
# the function `plot_PhytoMap` only can be run once. If you get an error plotting, please run the whole notebook :)
# 
# Press play in the map and check the whale movements over the hotspots of phytoplankton

# %%
# plot
Map = plot_PhytoMap(geodata, phytoplankton);
Map

# %%

save_PhytoMap(month, depth, hotspot_concentration)


