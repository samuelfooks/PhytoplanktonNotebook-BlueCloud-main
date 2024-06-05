## Whale Tracker using ARCO Blue Cloud Phytoplankton Data Products

This project is a modified version of the notebooks available at [gis4-wildlife/PhytoplanktonNotebook-BlueCloud](https://github.com/gis4-wildlife/PhytoplanktonNotebook-BlueCloud). It utilizes Analysis Ready Cloud Optimized(ARCO) versions of gridded chlorophyll data (a proxy of the total phytoplankton biomass) in combination with wildlife tracking data and other geospatial data to provide an interactive exploration of possible relationships between the movement of marine mammals and these parameters.  

### Overview

The project combines global open ocean three-dimensional (3D) gridded products of (1) chlorophyll a concentration (Chla), into a single Cloud Optimized Zarr dataset.  This is hosted on an S3 location. By accessing the data in this way, users can efficiently overlay each of the modeled phytoplankton months without the need to download large NetCDF datasets.

### Features

- Consolidated data into a single zarr dataset.
- Utilizes web map application for visual representation.
- Enables filtering of plankton rasters on a monthly basis.
- Simplifies overlaying whale movement and Marine Protected Areas.

### Platform

The project is hosted on the EDITO Datalab Platform, providing easy data sharing and on-demand cloud computing resources for data scientists and conservation enthusiasts.

### Usage

To use this project, you can install the required dependencies using one of the following methods:

#### Using environment.yml

mamba env update -n base -f environment.yml

#### using requirements.txt

pip install -r requirements.txt

### Acknowledgements

- **Brian R Vallejo:** [GitHub Profile](https://github.com/bryanvallejo16)
- **BlueCloud Phytoplankton NetCDFs:** [BlueCloud](https://blue-cloud.d4science.org/group/zoo-phytoplankton_eov/workspace)
- **MPA Dataset:** UNEP-WCMC and IUCN (2024), Protected Planet: The World Database on Protected Areas (WDPA) and World Database on Other Effective Area-based Conservation Measures (WD-OECM) [Online], May 2024, Cambridge, UK: UNEP-WCMC and IUCN. Available at: [Protected Planet](www.protectedplanet.net).



