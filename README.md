# Whale Tracker using ARCO Blue Cloud Phytoplankton Data Products

Altered version from the notebooks found on https://github.com/gis4-wildlife/PhytoplanktonNotebook-BlueCloud.

In this version we use the same products from BlueCloud PhytoPlanktonEOV

We combine all the outputs into a single zarr dataset, and put this onto an s3 location.  Loading this data onto the web map app allows us
to rapidly overlay each of modelled phytoplankton months without having to download large NetCDF datasets.

We then display the map in a similar application as used in gis4wildlife, but with all the plankton rasters onto a single page.  With the ability to filter the plankton rasters on a monthly basis. Instead of saving multiple html pages for each plankton raster.  We can then easier overlay the movement of the whales, and Marine Protected areas, in an easy to use visual representation. 

All this on the EDITO Datalab Platform, which allows for simple sharing of data, and on demand cloud computing resources for data scientists and conservation enthusiasts. 


Acknowledgements

Brian R Vallejo
https://github.com/bryanvallejo16

BlueCloud phytoplankton NetCDFs
https://blue-cloud.d4science.org/group/zoo-phytoplankton_eov/workspace	
Zoo-Phytoplankton_EOV >Phytoplankton_EOV >Chla_Product >Outputs >20210222 > 2018

MPA dataset:
UNEP-WCMC and IUCN (2024), Protected Planet: The World Database on Protected Areas (WDPA) and World Database on Other Effective Area-based Conservation Measures (WD-OECM) [Online], May 2024, Cambridge, UK: UNEP-WCMC and IUCN. Available at: www.protectedplanet.net.


