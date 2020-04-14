%% Data preparation for CAMELS data:
%   - Creates structure with CAMELS catchment data and catchment attributes
%   - Uses local paths (large CAMELS files)
%   - NA values in txt files are replaced with NaN values
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

% clc
% clear all
% close all

%% add directories for functions to path

%% load CAMELS data

% Before you can run the code you need to extract CAMELS data and store
% them in the correct folder. The following folders are required:
%   - camels_attributes_v2.0
%   - basin_timeseries_v1p2_modelOutput_daymet (used for time series)

disp('This function uses local paths. Change to local path where CAMELS data are stored.')

path_catchment_attributes = "C:\Users\sg16200\Local Documents\CAMELS_v2.0\camels_attributes_v2.0\"; % change to your local path
path_time_series = "C:\Users\sg16200\Local Documents\CAMELS_v2.0\basin_timeseries_v1p2_modelOutput_daymet\model_output_daymet\model_output\flow_timeseries\daymet\"; % change to your local path

if exist(path_catchment_attributes) == 7
    addpath(genpath(path_catchment_attributes)); % not actually needed since the files are just loaded
else
    error('Cannot find local path. You can download CAMELS from https://ral.ucar.edu/solutions/products/camels.')
end

%% load catchment attributes
% load data from files

% topography
% gauge_id;gauge_lat;gauge_lon;elev_mean;slope_mean;area_gages2;area_geospa_fabric
file_ID_topo = fopen(strcat(path_catchment_attributes,'camels_topo.txt'),'r');
file_topo = fread(file_ID_topo,'*char');
file_topo = strrep(file_topo','NA','NaN'); % replace NA with NaN
camels_topo_data = textscan(file_topo,'%f %f %f %f %f %f %f',...
    'Delimiter',';','headerlines',1);
fclose(file_ID_topo);

% climate
% gauge_id;p_mean;pet_mean;p_seasonality;frac_snow;aridity;high_prec_freq;high_prec_dur;high_prec_timing;low_prec_freq;low_prec_dur;low_prec_timing
file_ID_clim = fopen(strcat(path_catchment_attributes,'camels_clim.txt'),'r');
file_clim = fread(file_ID_clim,'*char');
file_clim = strrep(file_clim','NA','NaN'); % replace NA with NaN
camels_climate_data = textscan(file_clim,'%f %f %f %f %f %f %f %f %q %f %f %q',...
    'Delimiter',';','headerlines',1);
fclose(file_ID_clim);

% hydrology
% gauge_id;q_mean;runoff_ratio;slope_fdc;baseflow_index;stream_elas;q5;q95;high_q_freq;high_q_dur;low_q_freq;low_q_dur;zero_q_freq;hfd_mean
file_ID_hydro = fopen(strcat(path_catchment_attributes,'camels_hydro.txt'),'r');
file_hydro = fread(file_ID_hydro,'*char');
file_hydro = strrep(file_hydro','NA','NaN'); % replace NA with NaN
camels_hydro_data = textscan(file_hydro,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
    'Delimiter',';','headerlines',1);
fclose(file_ID_hydro);

% soils
% gauge_id;soil_depth_pelletier;soil_depth_statsgo;soil_porosity;soil_conductivity;max_water_content;sand_frac;silt_frac;clay_frac;water_frac;organic_frac;other_frac
file_ID_soil = fopen(strcat(path_catchment_attributes,'camels_soil.txt'),'r');
file_soil = fread(file_ID_soil,'*char');
file_soil = strrep(file_soil','NA','NaN'); % replace NA with NaN
camels_soil_data = textscan(file_soil,'%f %f %f %f %f %f %f %f %f %f %f %f',...
    'Delimiter',';','headerlines',1);
fclose(file_ID_soil);

% geology
% gauge_id;geol_1st_class;glim_1st_class_frac;geol_2nd_class;glim_2nd_class_frac;carbonate_rocks_frac;geol_porostiy;geol_permeability
file_ID_geol = fopen(strcat(path_catchment_attributes,'camels_geol.txt'),'r');
file_geol = fread(file_ID_geol,'*char');
file_geol = strrep(file_geol','NA','NaN'); % replace NA with NaN
camels_geol_data = textscan(file_geol,'%f %s %f %s %f %f %f %f',...
    'Delimiter',';','headerlines',1);
fclose(file_ID_geol);

% vegetation
% gauge_id;frac_forest;lai_max;lai_diff;gvf_max;gvf_diff;dom_land_cover_frac;dom_land_cover;root_depth_50;root_depth_99
file_ID_vege = fopen(strcat(path_catchment_attributes,'camels_vege.txt'),'r');
file_vege = fread(file_ID_vege,'*char');
file_vege = strrep(file_vege','NA','NaN'); % replace NA with NaN
camels_vege_data = textscan(file_vege,'%f %f %f %f %f %f %f %s %f %f',...
    'Delimiter',';','headerlines',1);
fclose(file_ID_vege);

%% load catchment attributes
% extract data from structures
gauge_id = camels_climate_data{:,1};

% topography
gauge_lat = camels_topo_data{:,2};
gauge_lon = camels_topo_data{:,3};
elev_mean = camels_topo_data{:,4};
slope_mean = camels_topo_data{:,5};
area_gages2 = camels_topo_data{:,6};
area_geospa_fabric = camels_topo_data{:,7};

% climate
p_mean = (camels_climate_data{:,2});
pet_mean = (camels_climate_data{:,3});
p_seasonality = (camels_climate_data{:,4});
frac_snow = (camels_climate_data{:,5});
aridity = (camels_climate_data{:,6});
high_prec_freq = (camels_climate_data{:,7});
high_prec_dur = (camels_climate_data{:,8});
high_prec_timing = (camels_climate_data{:,9});
low_prec_freq = (camels_climate_data{:,10});
low_prec_dur = (camels_climate_data{:,11});
low_prec_timing = (camels_climate_data{:,12});

% hydrology
q_mean = (camels_hydro_data{:,2});
runoff_ratio = (camels_hydro_data{:,3});
slope_fdc = (camels_hydro_data{:,4});
baseflow_index = (camels_hydro_data{:,5});
stream_elas = (camels_hydro_data{:,6});
q5 = (camels_hydro_data{:,7});
q95 = (camels_hydro_data{:,8});
high_q_freq = (camels_hydro_data{:,9});
high_q_dur = (camels_hydro_data{:,10});
low_q_freq = (camels_hydro_data{:,11});
low_q_dur = (camels_hydro_data{:,12});
zero_q_freq = (camels_hydro_data{:,13});
hfd_mean = (camels_hydro_data{:,14});

% soil
soil_depth_pelletier = (camels_soil_data{:,2});
soil_depth_statsgo = (camels_soil_data{:,3});
soil_porosity = (camels_soil_data{:,4});
soil_conductivity = (camels_soil_data{:,5});
max_water_content = (camels_soil_data{:,6});
sand_frac = (camels_soil_data{:,7});
silt_frac = (camels_soil_data{:,8});
clay_frac = (camels_soil_data{:,9});
water_frac = (camels_soil_data{:,10});
organic_frac = (camels_soil_data{:,11});
other_frac = (camels_soil_data{:,12});

% geology
geol_1st_class = (camels_geol_data{:,2});
glim_1st_class_frac = (camels_geol_data{:,3});
geol_2nd_class = (camels_geol_data{:,4});
glim_2nd_class_frac = (camels_geol_data{:,5});
carbonate_rocks_frac = (camels_geol_data{:,6});
geol_porosity = (camels_geol_data{:,7});
geol_permeability = (camels_geol_data{:,8});

% vegetation
frac_forest = (camels_vege_data{:,2});
lai_max = (camels_vege_data{:,3});
lai_diff = (camels_vege_data{:,4});
gvf_max = (camels_vege_data{:,5});
gvf_diff = (camels_vege_data{:,6});
dom_land_cover_frac = (camels_vege_data{:,7});
dom_land_cover = (camels_vege_data{:,8});
root_depth_50 = (camels_vege_data{:,9});
root_depth_99 = (camels_vege_data{:,10});

%% load time series
flow_perc_complete = NaN(length(gauge_id),1);
P = cell(length(gauge_id),1); % precipitation
PET = cell(length(gauge_id),1); % potential evapotranspiration
Q = cell(length(gauge_id),1); % streamflow
T = cell(length(gauge_id),1); % temperature

% loop over all catchments
for i = 1:length(gauge_id)
    
    if mod(i,10) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,length(gauge_id))
    end
    
    ID = gauge_id(i);
    [P, PET, Q, T] = loadCatchment_ID_CAMELS(ID,path_time_series);
    
    if sum(isnan(P(:,2))) > 0
        disp('NaN P');
    elseif sum(isnan(PET(:,2))) > 0
        disp('NaN PET');
    elseif sum(isnan(Q(:,2))) > 0
        disp('NaN Q');
    end
    
    flow_perc_complete(i) = 100*(1-sum(isnan(Q(:,2)))./length(Q(:,2)));
    
end

%% create structure

% topography
CAMELS_data.gauge_id = gauge_id;
CAMELS_data.gauge_lat = gauge_lat;
CAMELS_data.gauge_lon = gauge_lon;
CAMELS_data.elev_mean = elev_mean;
CAMELS_data.slope_mean = slope_mean;
CAMELS_data.area_gages2 = area_gages2;
CAMELS_data.area_geospa_fabric = area_geospa_fabric;

% climate
CAMELS_data.p_mean = p_mean;
CAMELS_data.pet_mean = pet_mean;
CAMELS_data.p_seasonality = p_seasonality;
CAMELS_data.frac_snow = frac_snow;
CAMELS_data.aridity = aridity;
CAMELS_data.high_prec_freq = high_prec_freq;
CAMELS_data.high_prec_dur = high_prec_dur;
CAMELS_data.high_prec_timing = high_prec_timing;
CAMELS_data.low_prec_freq = low_prec_freq;
CAMELS_data.low_prec_dur = low_prec_dur;
CAMELS_data.low_prec_timing = low_prec_timing;

% hydrology
CAMELS_data.q_mean = q_mean;
CAMELS_data.runoff_ratio = runoff_ratio;
CAMELS_data.slope_fdc = slope_fdc;
CAMELS_data.baseflow_index = baseflow_index;
CAMELS_data.stream_elas = stream_elas;
CAMELS_data.q5 = q5;
CAMELS_data.q95 = q95;
CAMELS_data.high_q_freq = high_q_freq;
CAMELS_data.high_q_dur = high_q_dur;
CAMELS_data.low_q_freq = low_q_freq;
CAMELS_data.low_q_dur = low_q_dur;
CAMELS_data.zero_q_freq = zero_q_freq;
CAMELS_data.hfd_mean = hfd_mean;

% soil
CAMELS_data.soil_depth_pelletier = soil_depth_pelletier;
CAMELS_data.soil_depth_statsgo = soil_depth_statsgo;
CAMELS_data.soil_porosity = soil_porosity;
CAMELS_data.soil_conductivity = soil_conductivity;
CAMELS_data.max_water_content = max_water_content;
CAMELS_data.sand_frac = sand_frac;
CAMELS_data.silt_frac = silt_frac;
CAMELS_data.clay_frac = clay_frac;
CAMELS_data.water_frac = water_frac;
CAMELS_data.organic_frac = organic_frac;
CAMELS_data.other_frac = other_frac;

% geology
CAMELS_data.geol_1st_class = geol_1st_class;
CAMELS_data.glim_1st_class_frac = glim_1st_class_frac;
CAMELS_data.geol_2nd_class = geol_2nd_class;
CAMELS_data.glim_2nd_class_frac = glim_2nd_class_frac;
CAMELS_data.carbonate_rocks_frac = carbonate_rocks_frac;
CAMELS_data.geol_porosity = geol_porosity;
CAMELS_data.geol_permeability = geol_permeability;

% vegetation
CAMELS_data.frac_forest = frac_forest;
CAMELS_data.lai_max = lai_max;
CAMELS_data.lai_diff = lai_diff;
CAMELS_data.gvf_max = gvf_max;
CAMELS_data.gvf_diff = gvf_diff;
CAMELS_data.dom_land_cover_frac = dom_land_cover_frac;
CAMELS_data.dom_land_cover = dom_land_cover;
CAMELS_data.root_depth_50 = root_depth_50;
CAMELS_data.root_depth_99 = root_depth_99;

% hydro-meteorological time series
CAMELS_data.flow_perc_complete = flow_perc_complete;
CAMELS_data.P = P;
CAMELS_data.PET = PET;
CAMELS_data.Q = Q;
CAMELS_data.T = T;

% save file to table
save('./CAMELS_Matlab/Data/CAMELS_data.mat','CAMELS_data')
