Phytoplankton EOV docs
Phytoplankton EOV products
 

The phytoplankton Essential Ocean Variables (EOV) demonstrator aims to provide a methodology to generate global open ocean three-dimensional (3D) gridded products of (1) chlorophyll a concentration (Chla), which is a proxy of the total phytoplankton biomass, and (2) Phytoplankton Functional Types (PFT), as a proxy for phytoplankton diversity, based on vertically-resolved in situ data of ocean physical properties (temperature and salinity) matched up with satellite products of ocean color and sea level anomaly.

 

The Machine Learning method
The methods have been developed following the method of Sauzede et al. (2016), which relies on machine learning, specifically on an artificial neural network (Multi-Layer Perceptron, MLP), and retrieves the vertical distribution of biogeochemical properties from merged ocean colour and hydrological data. The MLPs consist of several layers: one input layer, one output layer and one or several hidden layers. Each layer is composed of neurons, which are elementary transfer functions that provide outputs when inputs are applied.

Here, following the same philosophy as the method developed by Sauzede et al. (2016), two different MLP-based algorithms are developed for the independent retrieval of the Chla and of the PFT EOV products:

 

The first MLP retrieves the depth-resolved Chla product and is trained using in-situ depth-resolved measurements of Chla, temperature and salinity (T/S), from the global BioGeoChemical-Argo (BGC-Argo) observation network (coriolis Global Data Center), matched-up with global satellite-derived products. The MLP input layer is composed of three main components:
Surface satellite-based inputs from the Copernicus Marine Environment Monitoring Service (CMEMS) and GlobColour, such as the ocean colour remote sensing reflectance (Rrs) at five wavelengths, the photosynthetically available radiation (PAR), and the sea level anomaly (SLA),
Depth-resolved ocean physical properties such as components derived from a principal component analysis (PCA) of the T/S vertical profiles and the mixed layer depth ,
Time (day of the year transformed in cycles) and geographical coordinates of the ocean colour and hydrological data.
The retrieval of the depth-resolved PFT product relies on two distinct MLP-based algorithms. The first MLP is an upgraded version of the method developed by Sauzede et al. (2015) that includes hydrological information as input. This MLP is trained using a database comprising concurrent shipborne measurements of pigments determined by High Performance Liquid Chromatography (HPLC), chlorophyll fluorescence and T/S profiles. This MLP is applied to the BGC-Argo database, with vertical profiles of chlorophyll fluorescence and T/S used as inputs, in order to infer phytoplankton community composition, expressed in terms of depth-resolved Chla associated with three PFTs (pico-, nano- and microphytoplankton). This approach enables to enrich the global BGC-Argo database with the PFT information, which would not be available otherwise. Then, a second MLP is trained using the ‚PFT-enriched‚ BGC-Argo database, matched-up with satellite-derived products, in an analogous manner as described in 1.
 

Data sources
VARIABLES	DATA SOURCES	DATA ACCESS
Satellite-derived reflectance	OCEANCOLOUR_GLO_OPTICS_L3_REP_OBSERVATIONS_009_086	Blue Cloud
Satellite Sea Level Anomaly	SEALEVEL_GLO_PHY_L4_REP_OBSERVATIONS_008_047 product	Blue Cloud
Physical data: T, S, MLDGlobal ARMOR 3D products	CMEMS MULTIOBS_GLO_PHY_REP_015_002	Blue Cloud
BGC-Argo Float NetCDF files (S-files)	ftp.ifremer.fr/ifremer/argo/; http://www.argo.ucsd.edu	Blue Cloud
Satellite-derived Photosynthetically Available Radiation	ftp://ftp.hermes.acri.fr	Blue Cloud Vlab
High-performance liquid chromatography (HPLC) data	http://www.obs-vlfr.fr/proof/cruises.php	Blue Cloud Vlab
Bathymetry	https://www.gebco.net/data_and_products/gridded_bathymetry_data/	Blue Cloud Vlab
Data access= Blue Cloud, Data is accessible via the Blue Cloud Data Access service. It is also uploaded in the Vlab.

Data access= Blue Cloud Vlab, Data is not accessible via the Blue Cloud Data Access service. It has been uploaded in the Vlab.

 

 

Step by step guideline to use the service
All necessary files to generate monthly global 3D Chla product and the Phytoplankton functional types (PFT) for the year 2018 are located in Workspace/VRE Folders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV.

The input files are located in the Inputs folder in the Phytoplankton_EOV folder. The Inputs folder contains all the necessary input data to derive the Chla and PFT products for each month of the year 2018 (.nc files). Both products use the same input files.

Screenshot_folder
 

The Chla_Product and PFT_Product folders include three subfolders: Programs, Outputs, and Plots.

The Programs folder contains 2 Jupyter notebooks and two folders (Functions and Models).
The Functions folder contains all the necessary functions required to generate the 3D global products, and
The folder Models contains the trained MLP models and PCA models.
The Outputs folder contains the output global 3D products generated for each month of the year 2018 (.nc files).
The Plots folder contains the visualization of the output products as ‚ .png files (2D spatial plots for 36 depths).
 

 

How to run the Jupyter notebooks
To run the Jupyter notebooks users should make a copy of the folder Phytoplankton_EOV into their home directory (/home/jovyan/). To do this, launch the Jupyter hub and run the script in Workspace/VRE Folders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/copy_files_phyto.ipynb to automatically copy only the necessary files from the Phytoplankton_EOV folder. The Inputs can be called to the script directly from the VREFolders, and the Chla_Product and PFT_Product from the home directory.

Before executing these notebooks, verify the paths in the files: CREATE_MONTHLY_FIELDS_Loop_ZNORM.ipynb,Output_spatial_plots.ipynb, and Functions/SOCA_CHLA_ZNORM_2020.ipynb.

 

To generate the Chla product, open and run the two Jupyter notebooks available in the Phytoplankton_EOV/Chla_Product/Programs folder from the Jupyter Lab of the VRE, in the following order:
The first notebook: CREATE_MONTHLY_FIELDS_Loop_ZNORM.ipynb generates the global 3D Chla products in NetCDF format. For each month, the output is saved in the corresponding monthly folder, under Outputs.
The second notebook: Output_spatial_plots.ipynb is used to generate the visualization plots based on the output NetCDF files obtained from the first notebook. For each month, the plots are saved in the corresponding monthly folder under Plots.
To generate the PFT product, open and run the two Jupyter notebooks available in the Phytoplankton_EOV/PFT_Product/Programs folder from the Jupyter Lab of the VRE, in the following order:
The first notebook: CREATE_MONTHLY_FIELDS_PFT_ZNORM_N1.ipynb generates the global 3D PFT products (Micro-Chla, Nano-Chla, and Pico-Chla) in NetCDF format. For each month, the output is saved in the corresponding monthly folder, under Outputs.
The second notebook: Plots_output_spatial_monthly_PFT_2018.ipynb is used to generate the visualization plots based on the output NetCDF files obtained from the first notebook. For each month, the plots are saved in the corresponding monthly folder under Plots.
These notebooks generate monthly files for the 12 months of year 2018. The users can edit the variable MONTHS in the first notebook to create different products for other years. In the second notebook, the users can also plot different monthly png files by editing the variable MONTHS.

 

 

Using other data
Presently, this demonstrator generates global 3D Chla and PFT products for the year 2018. To generate the products for another year, users must have their input data in the same format as the data provided in the Inputs folder and edit the paths on the 2 main Jupyter notebooks provided in the Programs folder. The output NetCDF files and 2D spatial plots (.png format) will be generated in the corresponding monthly folders under Outputs and Plots.

The figures provided below illustrate the output product generated for the surface (0m depth) to 1000m depth for the months of January and July, which are respectively typical winter and summer months. The Chla product is in units of mg of chlorophyll a m-3; the PFT product provides the chlorophyll a concentration associated with the micro-, nano-, and picophytoplankton size classes in units of mg of chlorophyll a m-3.

 

January_2018_Chla

 

July_2018_Chla

 

January_2018_PFT

 

July_2018_PFT

 

 



 

Provide us with your feedback

 

Authors
Renosh Pannimpullath Remanan, Raphaelle Sauzede, Julia Uitz and Herve Claustre. Institut de la Mer de Villefranche, CNRS‚ Sorbonne Universite (France).

Maintainers: julia.uitz@imev-mer.fr and raphaelle.sauzede@imev-mer.fr

DOI (product): https://doi.org/10.48670/moi-00046

References
Sauzede, R., H. Claustre, C. Jamet, J. Uitz, J. Ras, A. Mignot, and F. d'Ortenzio (2015). Retrieving the vertical distribution of chlorophyll-a concentration and phytoplankton community composition from in situ fluorescence profiles: A method based on a neural network with potential for global-scale applications, J. Geophys. Res. Oceans, 120, 451-470, doi:10.1002/2014JC010355.
Sauzede, R., Claustre, H., Uitz, J., Jamet, C., Dall'Olmo, G., d'Ortenzio, F., Gentili, B., Poteau, A. and Schmechtig, C. (2016). A neural network-based method for merging ocean colour and Argo data to extend surface bio-optical properties to depth: Retrieval of the particulate backscattering coefficient. Journal of Geophysical Research: Oceans, 121(4), pp.2552-2571. http://doi.org/10.1002/2015JC011408
Cite as
The phytoplankton product developed in Blue Cloud has been used to update the Copernicus product with DOI http://doi.org/10.48670/moi-00046, that is published in the data catalogue from the E.U. Copernicus Marine Service Information.