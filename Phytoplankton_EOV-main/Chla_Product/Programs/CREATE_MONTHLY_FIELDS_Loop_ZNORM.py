# %% [markdown]
# # Notebook to create Monthly 3D SOCA-Chl

# %% [markdown]
# ## How to use the Chla demonstrator?
# 
# Chla_Product folder include three subfolders: Programs, Outputs, and Plots.
# 
# The Programs folder contains 2 Jupyter notebooks and two folders (Functions and Models).
# The Functions folder contains all the necessary functions required to generate the 3D global products and the folder Models contains the trained MLP models and PCA models.
# The Outputs folder contains the output global 3D products generated for each month of the year 2018 (“.nc” files.)
# The Plots folder contains the visualization of the output products as “.png” files (2D spatial plots for 36 depths).
# 
# The “Chla_Product” demonstrator can be executed in the Jupyter Lab of the VRE, by running the two Jupyter notebooks available in the “Phytoplankton_EOV/Chla_Product/Programs” folder, i.e. “CREATE_MONTHLY_FIELDS_Loop_ZNORM.ipynb” and “Output_spatial_plots.ipynb”. The first notebook generates the global 3D Chla products in NetCDF format. For each month, the output is saved in the corresponding monthly folder, under “Outputs”. The second notebook is used to generate the visualization plots based on the output NetCDF files obtained from the first notebook. For each month, the plots are saved in the corresponding monthly folder under “Plots”. Before executing these notebooks (“CREATE_MONTHLY_FIELDS_Loop_ZNORM.ipynb”, “Output_spatial_plots.ipynb”, and “Functions/SOCA_CHLA_ZNORM_2020.ipynb”), the paths should be checked and modified accordingly.
# 

# %% [markdown]
# ## Import libraries

# %%
#Import all functions that we need
import os
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset as NetCDFFile 
from datetime import datetime
from calendar import monthrange
import oceans.sw_extras as swe
import gsw
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore') # Ignore warning messages for printing

# %% [markdown]
# ## Import functions

# %% [markdown]
# ### Function to create the 3D Monthly NetCDF

# %%
import import_ipynb
from Functions.fonction_mask_black_sea import mask_black_sea
from Functions.fonction_mask_bathymetry import mask_bathy
from Functions.Function_find_Ze_depth import find_Ze
from Functions.SOCA_CHLA_ZNORM_2020 import INPUTS_SOCA_CHLA_2020, SOCA_CHLA_ZNORM_2020
from Functions.fonction_NetCDF_monthly_N import creation_NetCDF_3D_PRODUCT_CMEMS_N

# %% [markdown]
# ### Define paths for inputs and outputs
# 

# %%
# Define date and paths
# DATE_TODAY="20210222"
DATE_TODAY = '20211201' # DATE_TODAY is the Notebook running date; the output files will be saved under the folder Output with names DATE_TODAY

###############################################################################################
# # To run in home folder please uncomment below paths
path_TACMOB='/'.join(['/home/jovyan/Phytoplankton_EOV/Chla_Product/Outputs', DATE_TODAY])
path_data_sla = '/home/jovyan/Phytoplankton_EOV/Inputs/SLA'
path_data_rrs412 = '/home/jovyan/Phytoplankton_EOV/Inputs/RRS412'
path_data_rrs443 = '/home/jovyan/Phytoplankton_EOV/Inputs/RRS443'
path_data_rrs490 = '/home/jovyan/Phytoplankton_EOV/Inputs/RRS490'
path_data_rrs555 = '/home/jovyan/Phytoplankton_EOV/Inputs/RRS555'
path_data_rrs670 = '/home/jovyan/Phytoplankton_EOV/Inputs/RRS670'
path_data_par = '/home/jovyan/Phytoplankton_EOV/Inputs/PAR'
path_data_phy = '/home/jovyan/Phytoplankton_EOV/Inputs/ARMOR3D_N/'
path_data_chla = '/home/jovyan/Phytoplankton_EOV/Inputs/CHLA'
path_bathy_data="/home/jovyan/Phytoplankton_EOV/Inputs/BATHYMETRY/GEBCO_2014_6x6min_Global.nc"
###############################################################################################


###############################################################################################
# # To run in WorkSpace VRE folder please uncomment below paths
# path_TACMOB='/'.join(['/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Chla_Product/Outputs', DATE_TODAY])
# path_data_sla = '/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/SLA'
# path_data_rrs412 = '/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/RRS412'
# path_data_rrs443 = '/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/RRS443'
# path_data_rrs490 = '/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/RRS490'
# path_data_rrs555 = '/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/RRS555'
# path_data_rrs670 = '/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/RRS670'
# path_data_par = '/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/PAR'
# path_data_phy = '/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/ARMOR3D_N/'
# path_data_chla = '/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/CHLA'
# path_bathy_data="/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Inputs/BATHYMETRY/GEBCO_2014_6x6min_Global.nc"
###############################################################################################


# If the path where we will store the 3D NetCDF files is not a directory --> mkdir
if not os.path.isdir(path_TACMOB):
    os.mkdir(path_TACMOB)

# %% [markdown]
# ### **Important** - Input paths and FTP links to download monthly data
# 
# | Variable | Path_folder | Original data dowloaded from |
# | ---| --- | --- |
# | RRS412 | `path_data_rrs412` | ftp://my.cmems-du.eu/Core/OCEANCOLOUR_GLO_OPTICS_L4_REP_OBSERVATIONS_009_081/dataset-oc-glo-opt-multi-l4-rrs412_4km_monthly-rep-v02/ |
# | RRS443 | `path_data_rrs443` | ftp://my.cmems-du.eu/Core/OCEANCOLOUR_GLO_OPTICS_L4_REP_OBSERVATIONS_009_081/dataset-oc-glo-opt-multi-l4-rrs443_4km_monthly-rep-v02/ |
# | RRS490 | `path_data_rrs490` | ftp://my.cmems-du.eu/Core/OCEANCOLOUR_GLO_OPTICS_L4_REP_OBSERVATIONS_009_081/dataset-oc-glo-opt-multi-l4-rrs490_4km_monthly-rep-v02/ |
# | RRS555 | `path_data_rrs555` | ftp://my.cmems-du.eu/Core/OCEANCOLOUR_GLO_OPTICS_L4_REP_OBSERVATIONS_009_081/dataset-oc-glo-opt-multi-l4-rrs555_4km_monthly-rep-v02/ |
# | RRS670 | `path_data_rrs670` | ftp://my.cmems-du.eu/Core/OCEANCOLOUR_GLO_OPTICS_L4_REP_OBSERVATIONS_009_081/dataset-oc-glo-opt-multi-l4-rrs670_4km_monthly-rep-v02/ |
# | PAR |  `path_data_par` | ftp://ftp.hermes.acri.fr/GLOB/merged/month/ |
# | SLA |  `path_data_sla` | ftp://my.cmems-du.eu/Core/SEALEVEL_GLO_PHY_L4_REP_OBSERVATIONS_008_047/dataset-duacs-rep-global-merged-allsat-phy-l4-monthly/ |
# | Physical_ARMOR3D |  `path_data_phy` | ftp://nrt.cmems-du.eu/Core/MULTIOBS_GLO_PHY_TSUV_3D_MYNRT_015_012/dataset-armor-3d-rep-monthly/ |
# | CHLA |  `path_data_chla` | ftp://my.cmems-du.eu/Core/OCEANCOLOUR_GLO_OPTICS_L4_REP_OBSERVATIONS_009_081/dataset-oc-glo-opt-multi-l4-rrs490_4km_monthly-rep-v02/  |
# | BATHYMETRY |  `path_bathy_data` | https://download.gebco.net/ |
# 

# %% [markdown]
# ### Load bathymetry data

# %%
bathy_data=NetCDFFile(path_bathy_data)
Height=bathy_data.variables['Height'][:]
Lat_bathy=bathy_data.variables['lat'][:]
Lon_bathy=bathy_data.variables['lon'][:]
bathy_data.close()

# %% [markdown]
# ## Define a function to load a specific NetCDF file and return the matrix of parameter given as input with lon and lat

# %%
#Function to load data in NetCDF files
# Takes as input the path of the NetCDF file as well as the variable to return
def load_netcdf(path, var):
    
    #Open the NetCDF file and get the Chl/lon and latitude data
    nc=NetCDFFile(path)
    
    # Load data
    data=nc.variables[str(var)][:]
    
    # Get 2D matrix
    if len(data.shape)==4:
        data=data[0,0,:,:]
    if len(data.shape)==3:
        data=data[0,:,:]
        
    if var == 'sla':
        # SLA NetCDF have different variable names for longitude and latitude
        lon_data=nc.variables['longitude'][:]
        lon_data=list(lon_data)
        lat_data=nc.variables['latitude'][:]
        lat_data=list(lat_data)
        # Longitude range between 0 and 360 instead of -180 to 180 --> transformation for range -180-180
        for i in np.arange(len(lon_data)):
            #print(i)
            if lon_data[i]>=0 and lon_data[i]<180:
                lon_data[i]=lon_data[i]
            else:
                lon_data[i]=lon_data[i]-360

    else:
        lon_data=nc.variables['lon'][:]
        lon_data=list(lon_data)
        lat_data=nc.variables['lat'][:]
        lat_data=list(lat_data)
    
    
    #Then, close the NetCDF file
    nc.close()

    #Order the longitude vector and the Chl matrix
    sort_lon_data=np.argsort(lon_data)
    data=data[:,sort_lon_data]
    lon_data=np.array(lon_data)[sort_lon_data]

    #Idem for latitude
    sort_lat_data=np.argsort(lat_data)
    data=data[sort_lat_data,:]
    lat_data=np.array(lat_data)[sort_lat_data]

    return data, lon_data, lat_data

# %% [markdown]
# ## Define a function that returns the mean from different pixels

# %%
#Function to get the mean for different pixels in order to have the same resolution for satellite data as for physical data
def get_mean_pixel(VAR_MAT, lon_ix, lat_ix):
    #Take var alues for all iii and jjj --> mean value of these values is the matchup
    VAR=VAR_MAT[np.ix_(lat_ix, lon_ix)]
    #check if all values are masked: if yes put nan and if no compute the mean value
    if len(VAR[~VAR.mask])==0:
        MEAN=np.nan
    else:
        MEAN=np.nanmean(VAR[~VAR.mask])
    return MEAN

# %% [markdown]
# ### Compute ZNORM depth from surface Chla and MLD

# %%
def ZNORM_N(chl_sat,mld):
    Zeu = find_Ze(chl_sat)
    ZNORM = max(1.5*(Zeu),mld)
    return ZNORM

# %% [markdown]
# ### Create monthly .nc files

# %%
#vector depth contains all pressure of estimations (retrieval for the 3D cube, defined from ARMOR3D resolution)

CMEMS_pres_N=[0, 5, 10, 15, 20, 25, 30, 35, 40,45, 50, 55, 60, 65, 70, 80, 90, 100,
              125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, 600,
              700, 800, 900, 1000]
y = 2018
MONTHS=[1,2,3,4,5,6,7,8,9,10,11,12]
# By default this notebook generate monthly files for 12 months of year 2018; 
# If user want to generate the output for a special year, they can edit  MONTHS as follows 
# MONTHS=[1] 
# it generate product for January

str_months = ['January', 'February', "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
for m in MONTHS:
    print(y)
    print(m) 
    M = str(m).zfill(2)
    
    #Create the subdirectory of path_TACMOB where the NetCDF of 3D Chl will be stored
    
    dir_out_year = '/'.join([path_TACMOB,str(y)])
    if not os.path.isdir(dir_out_year):
        os.mkdir(dir_out_year)
    
    dir_nc_tac_mob=[path_TACMOB,str(y),str(m)]
    path_nc_tac_mob='/'.join(dir_nc_tac_mob)
    if not os.path.isdir(path_nc_tac_mob):
        os.mkdir(path_nc_tac_mob)
        
    dir_phy_data=[path_data_phy,str(y)]
    path_phy='/'.join(dir_phy_data)
    #go in this path to list every file with physical data (one file for one week of data)
    os.chdir(path_phy)
    nc_file_chemin_phy=os.listdir()
    
    #nc_file_chemin_phy=[sys.argv[2]]
#     f = ''.join(['dataset-armor-3d-rep-monthly_',str(y),str(M),'15T1200Z_P20190301T0000Z.nc'])
#     for f in nc_file_chemin_phy:
#     print(f)
    
    FFF = ''.join([str(y),str(M),'15'])
    f = ''.join([f1 for f1 in nc_file_chemin_phy if os.path.isfile(os.path.join(path_phy,f1)) and FFF in f1])
#     ''.join(f)
#     for f in nc_file_chemin_phy:
    print(f)
        
    #Beginning time to compute the time to create one NetCDF
    time1=datetime.now()
        
    #Get the date to have the corresponding ocean color data file
    date=f.split("_")[1][:8]
        
    #Get corresponding day of the year
    date_format_datetime =  datetime.strptime(date,'%Y%m%d')
    doy_f=float(date_format_datetime.strftime('%j'))
        
    path_NetCDF_phy=f
    #open the NetCDF with physical data and get the MLD/lon/lat matrix
    phy_data=NetCDFFile(path_NetCDF_phy)
    #Get mixed layer depth
    mlotst=phy_data.variables['mlotst'][:]
    #[0,:,:] 0 because the 1st dimension of the 3D matrix was 1 --> so transformation in 2D matrix
    mlotst=mlotst[0,:,:]
    # Get temperature
    to=phy_data.variables['to'][:]
    #[0,0:19,:,:] 0 because the 1st dimension of the 4D matrix was 1 --> so transformation in 3D matrix
    #0:19 to get the 19 depth for SOCA (from surface to 1000m depth)
    to=to[0,0:36,:,:]
    #Get salinity
    so=phy_data.variables['so'][:]
    #[0,0:19,:,:] 0 because the 1st dimension of the 4D matrix was 1 --> so transformation in 3D matrix
    #0:19 to get the 19 depth for SOCA (from surface to 1000m depth)
    so=so[0,0:36,:,:]
    #Get lon and lat
    lon_phy=phy_data.variables['longitude'][:]
    lon_phy=list(lon_phy)
    lat_phy=phy_data.variables['latitude'][:]
    lat_phy=list(lat_phy)
        
    #close the NetCDF physical file
    phy_data.close()      
        
    #transform longitude of physical file with range 0-360 with a new longitude with range -180 à 180 
    #to have the same range of longitude between phy and ocean color data --> for the matchup
    lon_phy2=lon_phy.copy()
    for i in np.arange(len(lon_phy2)):
            #print(i)
        if lon_phy[i]>=0 and lon_phy[i]<180:
            lon_phy2[i]=lon_phy2[i]
        else:
            lon_phy2[i]=lon_phy2[i]-360
                
    #Then order the longitude vector and the MLD/temp and sal matrix
    sort_lon_phy2=np.argsort(lon_phy2)
    mlotst=mlotst[:,sort_lon_phy2]
    to=to[:,:,sort_lon_phy2]
    so=so[:,:,sort_lon_phy2]
    lon_phy2=np.array(lon_phy2)[sort_lon_phy2]
        
    #Get the satellite data path and files names (with the good year and date)
    dir_data_rrs412='/'.join([path_data_rrs412,str(y)])
    num_days = monthrange(y, m)
    fin_day = num_days[1] # final day of the month       
        
    dir_data_rrs443='/'.join([path_data_rrs443,str(y)])
    dir_data_rrs490='/'.join([path_data_rrs490,str(y)])
    dir_data_rrs555='/'.join([path_data_rrs555,str(y)])
    dir_data_rrs670='/'.join([path_data_rrs670,str(y)])
    dir_data_chla='/'.join([path_data_chla,str(y)])
               
    dir_data_par='/'.join([path_data_par,str(y)])
    dir_data_sla='/'.join([path_data_sla,str(y)])
        
               
    path_NetCDF_rrs412=''.join([dir_data_rrs412,'/',str(y),str(M),'01_m_',str(y),str(M),str(fin_day),'-ACRI-L4-RRS412-AVW_MULTI_4KM-GLO-REP-v02.nc'])
    print(path_NetCDF_rrs412)    
        
    path_NetCDF_rrs443=''.join([dir_data_rrs443,'/',str(y),str(M),'01_m_',str(y),str(M),str(fin_day),'-ACRI-L4-RRS443-AVW_MULTI_4KM-GLO-REP-v02.nc'])
    print(path_NetCDF_rrs443)
        
    path_NetCDF_rrs490=''.join([dir_data_rrs490,'/',str(y),str(M),'01_m_',str(y),str(M),str(fin_day),'-ACRI-L4-RRS490-AVW_MULTI_4KM-GLO-REP-v02.nc'])
    print(path_NetCDF_rrs490)
        
    path_NetCDF_rrs555=''.join([dir_data_rrs555,'/',str(y),str(M),'01_m_',str(y),str(M),str(fin_day),'-ACRI-L4-RRS555-AVW_MULTI_4KM-GLO-REP-v02.nc'])
    print(path_NetCDF_rrs555)
        
    path_NetCDF_rrs670=''.join([dir_data_rrs670,'/',str(y),str(M),'01_m_',str(y),str(M),str(fin_day),'-ACRI-L4-RRS670-AVW_MULTI_4KM-GLO-REP-v02.nc'])
    print(path_NetCDF_rrs670)
        
    path_NetCDF_chla=''.join([dir_data_chla,'/',str(y),str(M),'01_m_',str(y),str(M),str(fin_day),'-ACRI-L4-CHL-MULTI_4KM-GLO-REP.nc'])
    print(path_NetCDF_chla)  
        
    path_NetCDF_par =''.join([dir_data_par,'/L3m_',str(y),str(M),'01-',str(y),str(M),str(fin_day),'__GLOB_4_AVW-MODVIR_PAR_MO_00.nc'])
    print(path_NetCDF_par)         
        
    path_NetCDF_sla=''.join([dir_data_sla,'/', 'dt_global_allsat_msla_h_y', str(y),'_m',str(M),'.nc'])
    print(path_NetCDF_sla)
    #####################################################################
    
    rrs412 = load_netcdf(path=path_NetCDF_rrs412, var="RRS412")[0]
    lon_oc = load_netcdf(path=path_NetCDF_rrs412, var="RRS412")[1]
    lat_oc = load_netcdf(path=path_NetCDF_rrs412, var="RRS412")[2]
    
    rrs443 = load_netcdf(path=path_NetCDF_rrs443, var="RRS443")[0]
       
    rrs490 = load_netcdf(path=path_NetCDF_rrs490, var="RRS490")[0]
        
    rrs555 = load_netcdf(path=path_NetCDF_rrs555, var="RRS555")[0]
        
    rrs670 = load_netcdf(path=path_NetCDF_rrs670, var="RRS670")[0]       
    
    chla = load_netcdf(path=path_NetCDF_chla, var="CHL")[0]        
    
    par = load_netcdf(path=path_NetCDF_par, var="PAR_mean")[0]
    
    
    sla = load_netcdf(path=path_NetCDF_sla, var="sla")[0]
    lon_sla = load_netcdf(path=path_NetCDF_sla, var="sla")[1]
    lat_sla = load_netcdf(path=path_NetCDF_sla, var="sla")[2]
        
    
       
        
    #Create VAR3D matrix which is the matrix in which we will have the cube retrieval of Chl from Uitz et al 2006
    SOCA_CHLA_3D=np.empty((len(CMEMS_pres_N),len(lat_phy),len(lon_phy2)))
    SOCA_CHLA_3D[:,:,:]=np.nan
    #Create CHL3D ERROR matrix which is the matrix in which we will have the cube error of retrieval of Chl from Uitz et al 2006
    SOCA_CHLA_3D_ERR=np.empty((len(CMEMS_pres_N),len(lat_phy),len(lon_phy2)))
    SOCA_CHLA_3D_ERR[:,:,:]=np.nan


    #Loop for each pixel (lon/lat phy)
    for i in np.arange(len(lon_phy2)):
        print(i)
        warnings.filterwarnings('ignore') # Ignore warning messages for printing
        #iii_oc the location of the pixel of ocean color that match with physical longitude
        iii_oc=np.logical_and(lon_oc >= lon_phy2[i]-0.125, lon_oc <= lon_phy2[i]+0.125)
        iii_oc=np.where(iii_oc)
        iii_oc=np.array(iii_oc)[0,:]
            
        #iii_sla is the location of the pixel of sla that match with physical longitude
        iii_sla=np.logical_and(lon_sla >= lon_phy2[i]-0.125, lon_sla <= lon_phy2[i]+0.125)
        iii_sla=np.where(iii_sla)
        iii_sla=np.array(iii_sla)[0,:]

        for j in np.arange(len(lat_phy)):
            #jjj_oc is the location of the pixel of ocean color that match with physical latitude
            jjj_oc=np.logical_and(lat_oc >= lat_phy[j]-0.125, lat_oc <= lat_phy[j]+0.125)
            jjj_oc=np.where(jjj_oc)
            jjj_oc=np.array(jjj_oc)[0,:]
                
            #jjj_sla is the location of the pixel of sla that match with physical latitude
            jjj_sla=np.logical_and(lat_sla >= lat_phy[j]-0.125, lat_sla <= lat_phy[j]+0.125)
            jjj_sla=np.where(jjj_sla)
            jjj_sla=np.array(jjj_sla)[0,:]
                
            # Apply the Black Sea mask:
            black_sea=mask_black_sea(lon=lon_phy2[i], lat=lat_phy[j])   
            # Apply the bathymetric mask <ith threshold 1500m depth:
            bathy=mask_bathy(bathy_mat=Height, lon_bathy=Lon_bathy, lat_bathy=Lat_bathy, lon=lon_phy2[i], lat=lat_phy[j], threshold=-1500)
            if not black_sea and bathy:  
#                    print("ok")
                #Get the MLD and Chl value for Chl vertical distribution
                mld_soca=mlotst[j,i]
#                     chl_sat_uitz=get_mean_pixel(VAR_MAT=chl, lon_ix=iii_oc, lat_ix=jjj_oc)

                # Get the temperature and salinity profiles
                temp_soca = to[:,j,i]
                sal_soca = so[:,j,i]
                
             
                # Get the RRS values and SLA and PAR
                rrs412_soca=get_mean_pixel(VAR_MAT=rrs412, lon_ix=iii_oc, lat_ix=jjj_oc)
                rrs443_soca=get_mean_pixel(VAR_MAT=rrs443, lon_ix=iii_oc, lat_ix=jjj_oc)
                rrs490_soca=get_mean_pixel(VAR_MAT=rrs490, lon_ix=iii_oc, lat_ix=jjj_oc)
                rrs555_soca=get_mean_pixel(VAR_MAT=rrs555, lon_ix=iii_oc, lat_ix=jjj_oc)
                rrs670_soca=get_mean_pixel(VAR_MAT=rrs670, lon_ix=iii_oc, lat_ix=jjj_oc)
                par_soca=get_mean_pixel(VAR_MAT=par, lon_ix=iii_oc, lat_ix=jjj_oc)
                sla_soca=get_mean_pixel(VAR_MAT=sla, lon_ix=iii_sla, lat_ix=jjj_sla)
                chla_matchup = get_mean_pixel(VAR_MAT=chla, lon_ix=iii_oc, lat_ix=jjj_oc)
                ZNORM=ZNORM_N(chla_matchup,mld_soca)
                


                    # Condition also on temp and sal because sometimes some depths are masked and not others so mld is not nan but some values in the T/S profiles are
#                     isinstance(mld_uitz,float) and
                if isinstance(mld_soca,float) and all([isinstance(T, float) for T in temp_soca]) and all([isinstance(S, float) for S in sal_soca]):
                    if not np.isnan(rrs412_soca) and not np.isnan(ZNORM) and not np.isnan(mld_soca) and not np.isnan(rrs443_soca) and not np.isnan(rrs490_soca) and not np.isnan(rrs555_soca) and not np.isnan(rrs670_soca) and not np.isnan(sla_soca) and not np.isnan(par_soca):
                        
                        ############ DERIVE SPICINESS ######################################
                        spici_soca = swe.spice(pd.DataFrame(sal_soca), pd.DataFrame(temp_soca), pd.DataFrame(CMEMS_pres_N))
                        spici_soca = np.squeeze(spici_soca)
                
                        ################## DERIVE BRUNT VAISALA FREQUENCY###################
                        #Compute Brünt-Väisälä Frequency squared (N²)
                        # Compute first the absolute salinity and conservative temperature
#                         Absolute_Sal = gsw.SA_from_SP(sal_soca,temp_soca,lon=lon_phy2[i], lat=lat_phy[j]) #help(gsw.SA_from_SP)
#                         Conservative_temp = gsw.CT_from_t(sal_soca, temp_soca, CMEMS_pres)
#                         # N squared and pres associated
#                         N2 = gsw.Nsquared(Absolute_Sal, Conservative_temp, CMEMS_pres, lat=lat_phy[j])[0]
#                         p_N2 = gsw.Nsquared(Absolute_Sal, Conservative_temp, CMEMS_pres, lat=lat_phy[j])[1]
#                         # # Then interpolate to have the good pressure 
#                         interp1=interp1d(p_N2,N2,kind='nearest',fill_value="extrapolate")
#                         brunt_soca=interp1(CMEMS_pres)
                        
                                                
                        Zeta_interp = np.linspace(0,1.5,50)
                        Zeta = CMEMS_pres_N/ZNORM
                        temp_soca_ZNORM = np.interp(Zeta_interp, Zeta,temp_soca)
                        sal_soca_ZNORM = np.interp(Zeta_interp, Zeta,sal_soca)
                        spici_soca_ZNORM = np.interp(Zeta_interp, Zeta,spici_soca)
#                         brunt_soca_ZNORM = np.interp(Zeta_interp, Zeta,brunt_soca)
                        
                                                
                        temp_soca_ZNORM = temp_soca_ZNORM.reshape(-1,1)
                        sal_soca_ZNORM = sal_soca_ZNORM.reshape(-1,1)
                        spici_soca_ZNORM = spici_soca_ZNORM.reshape(-1,1)
                                                
                
                        #SOCA2020 retrieval
                        soca_chla_2020_inputs=INPUTS_SOCA_CHLA_2020(RHO_WN_412=rrs412_soca*3.14, RHO_WN_443= rrs443_soca*3.14,
                                                                    RHO_WN_490=rrs490_soca*3.14, RHO_WN_555= rrs555_soca*3.14, RHO_WN_670= rrs670_soca*3.14,
                                                                    SLA= sla_soca, PAR=par_soca, MLD=mld_soca,sal= sal_soca_ZNORM, temp=temp_soca_ZNORM, spici=spici_soca_ZNORM,
                                                                    lon=lon_phy2[i], lat=lat_phy[j], doy=doy_f)

                        chla_soca = SOCA_CHLA_ZNORM_2020(soca_chla_2020_inputs,ZNORM,pres_new=CMEMS_pres_N)[0]
                        chla_soca_err = SOCA_CHLA_ZNORM_2020(soca_chla_2020_inputs,ZNORM,pres_new=CMEMS_pres_N)[1]

#                             POC_soca = POC_from_bbp(bbp_soca, bbp_soca_err)[0]
#                             POC_soca_err = POC_from_bbp(bbp_soca, bbp_soca_err)[1]

                            # Fill POC3D with SOCA2020 retrieval
    
                        SOCA_CHLA_3D[:,j,i]=chla_soca
                        SOCA_CHLA_3D_ERR[:,j,i]=chla_soca_err
#                             SOCA_CHLA_3D_N[:,j,i]=chla_soca
#                             SOCA_CHLA_3D_ERR_N[:,j,i]=chla_soca_err
            
    date1=f.split("_")[1][:8]
    DATE_PRODUCT=date1
    LONG=lon_phy2
    print(LONG.shape)
    LAT=np.array(lat_phy,)
    print(LAT.shape)
    DEPTH=np.array(CMEMS_pres_N)
    DEPTH=np.array(DEPTH,dtype="int16")
    print(DEPTH.shape)
    VAR4D=np.expand_dims(SOCA_CHLA_3D, axis=0)
    print(VAR4D.shape)
    VAR4D_ERR=np.expand_dims(SOCA_CHLA_3D_ERR, axis=0)
    print(VAR4D_ERR.shape)
#         VAR4D2=np.expand_dims(SOCA_CHLA_3D_N, axis=0)
#         VAR4D_ERR2=np.expand_dims(SOCA_CHLA_3D_ERR_N, axis=0)
        #Creation of the NetCDF
    creation_NetCDF_3D_PRODUCT_CMEMS_N(path=path_nc_tac_mob, date_today= DATE_TODAY, date_product=DATE_PRODUCT,depth_input=DEPTH,longitude_input=LONG, latitude_input=LAT,variable_4d_input1=VAR4D,variable_4d_error_input1=VAR4D_ERR) 
#                                          variable_4d_input2=VAR4D2,variable_4d_error_input2=VAR4D_ERR2)
        #End time to compute the time to create the NetCDF
    time2=datetime.now()
    #Compute the time
    difference =time2-time1
    print(difference)
#     text_file = open(''.join([str(path_TACMOB),str(f),'_time_computation.txt']), "w")
#     n = text_file.write(str(difference))
#     text_file.close()
########################
########################       
   

# %%



