# %% [markdown]
# # Function ZNORM SOCA2020-CHLA

# %% [markdown]
# ### Import libraries

# %%
import pathlib
import pickle
import numpy as np
import pandas as pd
import oceans.sw_extras.sw_extras as swe
from itertools import chain
import os
os.chdir(os.path.dirname(__file__))
# %% [markdown]
# ### Definition of paths

# %%

###############################################################################################
# # To run in home folder please uncomment below paths
MODEL_PATH = pathlib.Path("../Models/ZNORM/")
###############################################################################################


###############################################################################################
# # To run in WorkSpace VRE folder please uncomment below paths
# MODEL_PATH = pathlib.Path("/workspace/VREFolders/Zoo-Phytoplankton_EOV/Phytoplankton_EOV/Chla_Product/Programs/Models/ZNORM/")



# %% [markdown]
# ### Load the PCA model to transform temperature, salinity and spiciness profiles in principal components

# %%
# Lad the PCAs model
pkl_filename_temp_pca = ''.join([str(MODEL_PATH),str("/temp_pca_50_V1.pkl")])
pca_temp = pickle.load(open(pkl_filename_temp_pca, "rb"))

pkl_filename_sal_pca = ''.join([str(MODEL_PATH),str("/sal_pca_50_V1.pkl")])
pca_sal = pickle.load(open(pkl_filename_sal_pca, "rb"))

pkl_filename_spici_pca = ''.join([str(MODEL_PATH),str("/spici_pca_50_V1.pkl")])
pca_spici = pickle.load(open(pkl_filename_spici_pca, "rb"))

pkl_filename_chla_pca = ''.join([str(MODEL_PATH),str("/chla_pca_50_V1.pkl")])
pca_chla = pickle.load(open(pkl_filename_chla_pca, "rb"))



# %% [markdown]
# ### Load the x-scaler and y-scaler models

# %%
# Load the x-scaler from the file 
pkl_filename_x_scaler = ''.join([str(MODEL_PATH),str("/x_scaler_50_V1.pkl")]) 
x_scaler = pickle.load(open(pkl_filename_x_scaler, 'rb')) 
# Load the x-scaler from the file 
pkl_filename_y_scaler = ''.join([str(MODEL_PATH),str("/y_scaler_50_V1.pkl")])
y_scaler = pickle.load(open(pkl_filename_y_scaler, 'rb'))  

# %% [markdown]
# ### Function SOCA2020 that takes as inputs the vector of Rrs data, the SLA value, the PAR value, the vectors of salinity and temperature 

# %%
def INPUTS_SOCA_CHLA_2020(RHO_WN_412, RHO_WN_443, RHO_WN_490, RHO_WN_555, RHO_WN_670, SLA, PAR, MLD, sal,temp,spici,lon, lat, doy):
#                           sin_doy, cos_doy, x_cart, y_cart, z_cart, sal, temp, spici):
    
    # Transform in pandas dataframe
#     temp = pd.DataFrame(temp)
#     sal = pd.DataFrame(sal)
#     pres = pd.DataFrame(pres)
    INPUTS1 = [SLA, PAR, RHO_WN_412, RHO_WN_443, RHO_WN_490, RHO_WN_555, RHO_WN_670, MLD]
#               sin_doy,cos_doy,x_cart,y_cart,z_cart]
    
    # Compute spiciness from temperature and salinity
#     spiciness = swe.spice(sal, temp, pres)

    # Apply PCA on temperature, salinity and spiciness profiles to get the principal components
    PrincipalComponentsTemp_Test = pca_temp.transform(np.transpose(temp))[:,:4]
    principalComponentsSal_Test = pca_sal.transform(np.transpose(sal))[:,:3]
    principalComponentsSpici_Test = pca_spici.transform(np.transpose(spici))[:,:4]
    
    # Transform day of the year in sin and cos of the radians
    doy_radians = (doy*np.pi)/182.5
    DOY = [np.sin(doy_radians), np.cos(doy_radians)]
    
    # Transform the locationlon/lat in Cartesian coordinates x/y/z
    lat_radians = (lat*np.pi)/90
    lon_radians = (lon*np.pi)/180
    X_Y_Z = [np.cos(lat_radians) * np.cos(lon_radians), np.cos(lat_radians) * np.sin(lon_radians), np.sin(lat_radians)]
         
    INPUTS = INPUTS1 + DOY + X_Y_Z + np.concatenate(PrincipalComponentsTemp_Test).ravel().tolist() + np.concatenate(principalComponentsSal_Test).ravel().tolist() + np.concatenate(principalComponentsSpici_Test).ravel().tolist()
    
    input_features = ['sla', 'PAR','RHO_WN_412', 
                  'RHO_WN_443','RHO_WN_490','RHO_WN_555', 'RHO_WN_670', 'MLD', 
                  'sin_doy', 'cos_doy','x_cart', 'y_cart', 'z_cart', 
                  'temp0', 'temp1', 'temp2', 'temp3', 'sal0', 'sal1', 'sal2', 'spici0', 'spici1', 'spici2', 'spici3']

#     input_features = ['sla', 'PAR', 
#                       'RHO_WN_412', 'RHO_WN_443','RHO_WN_490', 'RHO_WN_555', 'RHO_WN_670', 
#                       'MLD', 'doy_sin', 'doy_cos','x_cart', 'y_cart', 'z_cart', 
#                       'temp1', 'temp2', 'temp3', 'temp4', 'temp5', 'temp6', 'temp7',
#                       'sal1', 'sal2', 'sal3', 'spici1', 'spici2', 'spici3', 'spici4']


    INPUTS = pd.DataFrame(INPUTS)
    INPUTS = INPUTS.transpose()    
    INPUTS.columns = input_features
    INPUTS_scaled = pd.DataFrame(x_scaler.transform(INPUTS), columns=INPUTS.columns)
    
    return INPUTS_scaled

# %%
def SOCA_CHLA_ZNORM_2020(INPUTS_SOCA_CHLA_2020, ZNORM, pres_new):
    
    n_sample = INPUTS_SOCA_CHLA_2020.shape[0]
    YPRED_TOTAL = pd.DataFrame()
    nc1 = [65, 74, 76, 81, 86, 87, 92, 96, 98, 99]
    nc2 = [94, 73, 84, 83, 59, 84, 45, 99, 46, 72]
#     nc1 = [41, 44, 47, 48, 49, 49, 49, 50, 50, 50]
#     nc2 = [39, 34, 46, 48, 45, 48, 49, 46, 49, 50]
    Zeta_interp = np.linspace(0,1.5,50)
    depth_interp = Zeta_interp*ZNORM   
    
    
    for x in np.arange(len(nc1)):
        pkl_filename = ''.join([str(MODEL_PATH),str("/pickle_model_chllogpca50_V1_"),str(x),str(".pkl")])
        MLP = pickle.load(open(pkl_filename, 'rb'))
 
        # Apply MLP to the test dataset
        ypred_scaled = MLP.predict(INPUTS_SOCA_CHLA_2020)
        
        # Detransform the outputs:
        ypred = pd.DataFrame(
            y_scaler.inverse_transform(ypred_scaled),
            columns=np.arange(25),
            index=np.arange(n_sample)
        )

        # make vectors
        CHLA_PRED = pca_chla.inverse_transform(ypred)[:,:50]
        CHLA_PRED = pd.DataFrame(CHLA_PRED)
        CHLA_PRED_vec = CHLA_PRED.values.ravel()  
        CHLA_PRED_vec = pow(10,CHLA_PRED_vec)
        
        CHLA_PRED_vec_N = np.interp(pres_new, depth_interp,CHLA_PRED_vec) # made into 226 depths
        CHLA_PRED_vec_N = np.array(CHLA_PRED_vec_N)
        CHLA_PRED_vec_N[(pres_new > max(depth_interp))]=np.nan # masked data below 1.5 ZNORM depths              


        YPRED_TOTAL = pd.concat([YPRED_TOTAL, pd.DataFrame(CHLA_PRED_vec_N)], axis=1)
    
    soca_chla=np.nanmean(YPRED_TOTAL, axis=1)
    soca_chla_err=np.nanstd(YPRED_TOTAL, axis=1)
    
    return soca_chla, soca_chla_err

# %%



