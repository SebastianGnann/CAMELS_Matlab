%% saveCAMELSdata - saves various CAMELS datasets as struct files
%
%   This script loads various CAMELS datasets and saves them as Matlab 
%   struct files for easy use in Matlab. The timeseries are loaded only for 
%   the period in which all data are available (typically, streamflow data
%   are the shortest time series and therefore all other time series are
%   shortened to match the streamflow time series). Note that this script
%   can be slow and requires sufficient RAM.
%   
%   Links to the different CAMELS datasets are given below and instructions
%   where to place the datasets are given below. 
%
%   References
%
%   Newman, A.J., Clark, M.P., Sampson, K., Wood, A., Hay, L.E., Bock, A., 
%   Viger, R.J., Blodgett, D., Brekke, L., Arnold, J.R. and Hopson, T., 
%   2015. Development of a large-sample watershed-scale hydrometeorological 
%   data set for the contiguous USA: data set characteristics and 
%   assessment of regional variability in hydrologic model performance. 
%   Hydrology and Earth System Sciences, 19(1), p.209.
%   Addor, N., Newman, A.J., Mizukami, N. and Clark, M.P., 2017. The CAMELS
%   data set: catchment attributes and meteorology for large-sample
%   studies. Hydrology and Earth System Sciences (HESS), 21(10), 
%   pp.5293-5313.
%
%   Alvarez-Garreton, C., Mendoza, P.A., Boisier, J.P., Addor, N., 
%   Galleguillos, M., Zambrano-Bigiarini, M., Lara, A., Puelma, C., Cortes, 
%   G., Garreaud, R. and McPhee, J., 2018. The CAMELS-CL dataset: catchment 
%   attributes and meteorology for large sample studies–Chile dataset. 
%   Hydrology and Earth System Sciences, 22(11), pp.5817-5846.
%
%   Coxon, G., Addor, N., Bloomfield, J.P., Freer, J., Fry, M., Hannaford, 
%   J., Howden, N.J., Lane, R., Lewis, M., Robinson, E.L. and Wagener, T., 
%   2020. CAMELS-GB: hydrometeorological time series and landscape 
%   attributes for 671 catchments in Great Britain. Earth System Science 
%   Data, 12(4), pp.2459-2483.
%
%   Chagas, V.B., Chaffe, P.L., Addor, N., Fan, F.M., Fleischmann, A.S., 
%   Paiva, R.C. and Siqueira, V.A., 2020. CAMELS-BR: hydrometeorological 
%   time series and landscape attributes for 897 catchments in Brazil. 
%   Earth System Science Data, 12(3), pp.2075-2096.
%
%   Fowler, K.J., Acharya, S.C., Addor, N., Chou, C. and Peel, M.C., 2021. 
%   CAMELS-AUS: Hydrometeorological time series and landscape attributes 
%   for 222 catchments in Australia. Earth System Science Data Discussions, 
%   pp.1-30.
%
%   Copyright (C) 2021
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

close all
% clear all
clc

%% Data location and directories
% The code assumes that all the data are stored in the directory that we
% are in when we run the code. That is, when we run the code while being
% in the directory that contains CAMELS_Matlab (this repository), then the 
% CAMELS data files should also be in this directory. For example, if we 
% work in the main Matlab directory (here MATLAB), our folder structure 
% might look as follows:
%
% MATLAB/... (this is the directory we are in)
% MATLAB/CAMELS_Matlab (this repository contains the code we want to run)
% MATLAB/CAMELS_US (this folder contains CAMELS US data)
%                 
% If we want to use a different folder structure, we have to adjust the
% paths. 
% We need to add the CAMELS_Matlab repository to the Matlab path, so that
% we can work with relative paths. 
mydir = 'CAMELS_Matlab';
addpath(genpath(mydir));

%% CAMELS US
% First, we need to download and extract the CAMELS US data from:
% https://ral.ucar.edu/solutions/products/camels
% From the link above we can download six different folders. We need 3 
% of them to run this script:
% basin_timeseries_v1p2_modelOutput_daymet (CAMELS time series meteorology, observed flow, meta data (.zip))
% basin_timeseries_v1p2_metForcing_obsFlow (CAMELS time series Daymet forced model output (.zip))
% camels_attributes_v2.0 (CAMELS Attributes (.zip))
% The data should be stored in a folder named CAMELS_US located in the same
% directory as CAMELS_Matlab.

% Loading and saving the data might take a few minutes.
save_struct = true;
saveCAMELSstruct_US(save_struct); 

% Alternatively, we could also just load the struct file without saving it.
% save_struct = false;
% CAMELS_US_data = saveCAMELSstruct_US(save_struct); 

%% CAMELS GB 
% First, we need to download and extract the CAMELS GB data from:
% https://catalogue.ceh.ac.uk/documents/8344e4f3-d2ea-44f5-8afa-86d2987543a9
% From the link above we can download the CAMELS GB data all at once. 
% The data should be stored in a folder named CAMELS_GB located in the same
% directory as CAMELS_Matlab.

% Loading and saving the data might take a few minutes.
save_struct = true;
saveCAMELSstruct_GB(save_struct);

% Alternatively, we could also just load the struct file without saving it.
% save_struct = false;
% CAMELS_GB_data = saveCAMELSstruct_GB(save_struct); 

%% CAMELS CL 
% First, we need to download and extract the CAMELS CL data from:
% https://doi.pangaea.de/10.1594/PANGAEA.894885?format=html#download
% From the link above we can download 15 different folders. We need 7 of
% them to run this script (each contain one txt file):
% CAMELS_CL/1_CAMELScl_attributes
% CAMELS_CL/3_CAMELScl_streamflow_mm
% CAMELS_CL/5_CAMELScl_precip_chirps
% CAMELS_CL/6_CAMELScl_precip_mswep
% CAMELS_CL/10_CAMELScl_tmean_cr2met
% CAMELS_CL/11_CAMELScl_pet_8d_modis
% CAMELS_CL/12_CAMELScl_pet_hargreaves
% The data should be stored in a folder named CAMELS_CL located in the same
% directory as CAMELS_Matlab.

% Loading and saving the data might take a few minutes.
save_struct = true;
saveCAMELSstruct_CL(save_struct);

% Alternatively, we could also just load the struct file without saving it.
% save_struct = false;
% CAMELS_CL_data = saveCAMELSstruct_CL(save_struct); 

%% CAMELS BR 
% First, we need to download and extract the CAMELS BR data from:
% https://zenodo.org/record/3964745
% From the link above we can download 15 different folders. We need 5 of
% them to run this script (each contain multiple files):
% CAMELS_BR/1_CAMELS_BR_attributes
% CAMELS_BR/3_CAMELS_BR_streamflow_mm
% CAMELS_BR/4_CAMELS_BR_precipitation_chirps
% CAMELS_BR/9_CAMELS_BR_potential_evapotransp_gleam
% CAMELS_BR/11_CAMELS_BR_temperature_mean_cpc
% The data should be stored in a folder named CAMELS_BR located in the same
% directory as CAMELS_Matlab.

% Loading and saving the data might take a few minutes.
save_struct = true;
saveCAMELSstruct_BR(save_struct);

% Alternatively, we could also just load the struct file without saving it.
% save_struct = false;
% CAMELS_BR_data = saveCAMELSstruct_BR(save_struct); 


%% CAMELS AUS 
% First, we need to download and extract the CAMELS AUS data from:
% https://doi.pangaea.de/10.1594/PANGAEA.921850?format=html#download
% From the link above we can download 5 different folders and a few more 
% files. We need 3 folders/files to run this script:
% CAMELS_AUS/03_streamflow (contains streamflow time series)
% CAMELS_AUS/05_hydrometeorology (contains meteorological time series)
% CAMELS_AUS_Attributes-Indices_MasterTable (contains catchment attributes)
% The data should be stored in a folder named CAMELS_AUS located in the 
% same directory as CAMELS_Matlab.

% Loading and saving the data might take a few minutes.
save_struct = true;
saveCAMELSstruct_AUS(save_struct);

% Alternatively, we could also just load the struct file without saving it.
% save_struct = false;
% CAMELS_AUS_data = saveCAMELSstruct_AUS(save_struct); 

