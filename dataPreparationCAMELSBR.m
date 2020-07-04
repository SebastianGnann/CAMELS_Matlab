%% Data preparation for CAMELS BR data:
%   - Creates structure with CAMELS catchment data and catchment attributes
%   - Uses local paths (large CAMELS files)
%   - NA values in txt files are replaced with NaN values
%   - Loads CHIRPS P and GLEAM PET data
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
% them in the correct folder. The following folders/files are required:
%   \CAMELS_BR
%       \1_CAMELS_BR_attributes
%   	\3_CAMELS_BR_streamflow_mm
%   	\4_CAMELS_BR_precipitation_chirps
%   	\9_CAMELS_BR_potential_evapotransp_gleam
%   	\11_CAMELS_BR_temperature_mean_cpc

disp('This function uses local paths. Change to local path where CAMELS data are stored.')

path_catchment_attributes = "C:\Users\sg16200\Local Documents\CAMELS_BR\1_CAMELS_BR_attributes\"; % change to your local path
path_time_series = "C:\Users\sg16200\Local Documents\CAMELS_BR\"; % change to your local path

if ~exist(path_time_series) == 7
    error('Cannot find local path. You can download CAMELS BR from https://zenodo.org/record/3709338.')
end

%% load catchment attributes
% load data from files

% location
% gauge_id gauge_name gauge_region gauge_lat gauge_lon
file_ID_loca = fopen(strcat(path_catchment_attributes,'camels_br_location.txt'),'r');
file_loca = fread(file_ID_loca,'*char');
camels_BR_loca_data = textscan(file_loca,'%f %s %s %f %f',...
    'Delimiter',' ','headerlines',1);
fclose(file_ID_loca);

% topography
% gauge_id elev_gauge elev_mean slope_mean area
file_ID_topo = fopen(strcat(path_catchment_attributes,'camels_br_topography.txt'),'r');
file_topo = fread(file_ID_topo,'*char');
camels_BR_topo_data = textscan(file_topo,'%f %f %f %f %f',...
    'Delimiter',' ','headerlines',1);
fclose(file_ID_topo);

% climate
% gauge_id p_mean pet_mean aridity p_seasonality frac_snow high_prec_freq high_prec_dur high_prec_timing low_prec_freq low_prec_dur low_prec_timing
file_ID_clim = fopen(strcat(path_catchment_attributes,'camels_br_climate.txt'),'r');
file_clim = fread(file_ID_clim,'*char');
camels_BR_climate_data = textscan(file_clim,'%f %f %f %f %f %f %f %f %s %f %f %s',...
    'Delimiter',' ','headerlines',1);
fclose(file_ID_clim);

% hydrology
% gauge_id q_mean runoff_ratio stream_elas slope_fdc baseflow_index hfd_mean Q5 Q95 high_q_freq high_q_dur low_q_freq low_q_dur zero_q_freq
file_ID_hydro = fopen(strcat(path_catchment_attributes,'camels_br_hydrology.txt'),'r');
file_hydro = fread(file_ID_hydro,'*char');
camels_BR_hydro_data = textscan(file_hydro,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
    'Delimiter',' ','headerlines',1);
fclose(file_ID_hydro);

% soils
% gauge_id sand_frac silt_frac clay_frac org_carbon_content depth_bedrock water_table_depth
file_ID_soil = fopen(strcat(path_catchment_attributes,'camels_br_soil.txt'),'r');
file_soil = fread(file_ID_soil,'*char');
camels_BR_soil_data = textscan(file_soil,'%f %f %f %f %f %f %f',...
    'Delimiter',' ','headerlines',1);
fclose(file_ID_soil);

% geology
% gauge_id geol_class_1st geol_class_1st_frac geol_class_2nd geol_class_2nd_frac carb_rocks_frac geol_porosity geol_permeability
file_ID_geol = fopen(strcat(path_catchment_attributes,'camels_br_geology.txt'),'r');
file_geol = fread(file_ID_geol,'*char');
camels_BR_geol_data = textscan(file_geol,'%f %s %f %s %f %f %f %f',...
    'Delimiter',' ','headerlines',1);
fclose(file_ID_geol);

% land cover
% gauge_id crop_frac crop_mosaic_frac forest_frac shrub_frac grass_frac bare_frac imperv_frac wet_frac snow_frac dom_land_cover dom_land_cover_frac
file_ID_land = fopen(strcat(path_catchment_attributes,'camels_br_land_cover.txt'),'r');
file_land = fread(file_ID_land,'*char');
camels_BR_land_data = textscan(file_land,'%f %f %f %f %f %f %f %f %f %f %s %f',...
    'Delimiter',' ','headerlines',1);
fclose(file_ID_land);

% human intervention
% gauge_id consumptive_use_mm consumptive_use vol_reservoirs degree_of_regulation
file_ID_huma = fopen(strcat(path_catchment_attributes,'camels_br_human_intervention.txt'),'r');
file_huma = fread(file_ID_huma,'*char');
camels_BR_huma_data = textscan(file_huma,'%f %f %f %f %f',...
    'Delimiter',' ','headerlines',1);
fclose(file_ID_huma);

%quality check 
% gauge_id q_quality_control_frac q_stream_stage_frac
file_ID_qual = fopen(strcat(path_catchment_attributes,'camels_br_quality_check.txt'),'r');
file_qual = fread(file_ID_qual,'*char');
camels_BR_qual_data = textscan(file_qual,'%f %f %f',...
    'Delimiter',' ','headerlines',1);
fclose(file_ID_qual);

%% load catchment attributes
% extract data from structures
gauge_id = camels_BR_climate_data{:,1};

% location
gauge_name = camels_BR_loca_data{:,2};
gauge_region = camels_BR_loca_data{:,3};
gauge_lat = camels_BR_loca_data{:,4};
gauge_lon = camels_BR_loca_data{:,5};

% topography
elev_gauge = camels_BR_topo_data{:,2};
elev_mean = camels_BR_topo_data{:,3};
slope_mean = camels_BR_topo_data{:,4};
area = camels_BR_topo_data{:,5};

% climate
p_mean = (camels_BR_climate_data{:,2});
pet_mean = (camels_BR_climate_data{:,3});
aridity = (camels_BR_climate_data{:,4});
p_seasonality = (camels_BR_climate_data{:,5});
frac_snow = (camels_BR_climate_data{:,6});
high_prec_freq = (camels_BR_climate_data{:,7});
high_prec_dur = (camels_BR_climate_data{:,8});
high_prec_timing = (camels_BR_climate_data{:,9});
low_prec_freq = (camels_BR_climate_data{:,10});
low_prec_dur = (camels_BR_climate_data{:,11});
low_prec_timing = (camels_BR_climate_data{:,12});

% hydrology 
q_mean = (camels_BR_hydro_data{:,2});
runoff_ratio = (camels_BR_hydro_data{:,3});
stream_elas = (camels_BR_hydro_data{:,4});
slope_fdc = (camels_BR_hydro_data{:,5});
baseflow_index = (camels_BR_hydro_data{:,6});
hfd_mean = (camels_BR_hydro_data{:,7});
q5 = (camels_BR_hydro_data{:,8});
q95 = (camels_BR_hydro_data{:,9});
high_q_freq = (camels_BR_hydro_data{:,10});
high_q_dur = (camels_BR_hydro_data{:,11});
low_q_freq = (camels_BR_hydro_data{:,12});
low_q_dur = (camels_BR_hydro_data{:,13});
zero_q_freq = (camels_BR_hydro_data{:,14});

% soil
sand_frac = (camels_BR_soil_data{:,2});
silt_frac = (camels_BR_soil_data{:,3});
clay_frac = (camels_BR_soil_data{:,4});
org_carbon_content = (camels_BR_soil_data{:,5});
depth_bedrock = (camels_BR_soil_data{:,6});
water_table_depth = (camels_BR_soil_data{:,7});

% geology
geol_1st_class = (camels_BR_geol_data{:,2});
glim_1st_class_frac = (camels_BR_geol_data{:,3});
geol_2nd_class = (camels_BR_geol_data{:,4});
glim_2nd_class_frac = (camels_BR_geol_data{:,5});
carbonate_rocks_frac = (camels_BR_geol_data{:,6});
geol_porosity = (camels_BR_geol_data{:,7});
geol_permeability = (camels_BR_geol_data{:,8});

% land cover 
crop_frac = (camels_BR_land_data{:,2});
crop_mosaic_frac = (camels_BR_land_data{:,3});
forest_frac = (camels_BR_land_data{:,4});
shrub_frac = (camels_BR_land_data{:,5});
grass_frac = (camels_BR_land_data{:,6});
bare_frac = (camels_BR_land_data{:,7});
imperv_frac = (camels_BR_land_data{:,8});
wet_frac = (camels_BR_land_data{:,9});
snow_frac = (camels_BR_land_data{:,10});
dom_land_cover = (camels_BR_land_data{:,11});
dom_land_cover_frac = (camels_BR_land_data{:,12});

% human intervention
consumptive_use_mm = (camels_BR_huma_data{:,2});
consumptive_use = (camels_BR_huma_data{:,3});
vol_reservoirs = (camels_BR_huma_data{:,4});
degree_of_regulation = (camels_BR_huma_data{:,5});

% quality check
q_quality_control_frac = (camels_BR_qual_data{:,2});
q_stream_stage_frac = (camels_BR_qual_data{:,3});

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
    [P{i}, PET{i}, Q{i}, T{i}] = loadCatchmentCAMELSBR(ID,path_time_series);
    flow_perc_complete(i) = 100*(1-sum(isnan(Q{i}(:,2)))./length(Q{i}(:,2)));
    
end

%% create structure

% location
% TODO: extract 897 catchments
CAMELS_BR_data.gauge_id = gauge_id;
CAMELS_BR_data.gauge_name = gauge_name;
CAMELS_BR_data.gauge_region = gauge_region;
CAMELS_BR_data.gauge_lat = gauge_lat;
CAMELS_BR_data.gauge_lon = gauge_lon;

% topography
CAMELS_BR_data.elev_gauge = elev_gauge;
CAMELS_BR_data.elev_mean = elev_mean;
CAMELS_BR_data.slope_mean = slope_mean;
CAMELS_BR_data.area = area;

% climate
CAMELS_BR_data.p_mean = p_mean;
CAMELS_BR_data.pet_mean = pet_mean;
CAMELS_BR_data.aridity = aridity;
CAMELS_BR_data.p_seasonality = p_seasonality;
CAMELS_BR_data.frac_snow = frac_snow;
CAMELS_BR_data.high_prec_freq = high_prec_freq;
CAMELS_BR_data.high_prec_dur = high_prec_dur;
CAMELS_BR_data.high_prec_timing = high_prec_timing;
CAMELS_BR_data.low_prec_freq = low_prec_freq;
CAMELS_BR_data.low_prec_dur = low_prec_dur;
CAMELS_BR_data.low_prec_timing = low_prec_timing;

% hydrology 
CAMELS_BR_data.q_mean = q_mean;
CAMELS_BR_data.runoff_ratio = runoff_ratio;
CAMELS_BR_data.stream_elas = stream_elas;
CAMELS_BR_data.slope_fdc = slope_fdc;
CAMELS_BR_data.baseflow_index = baseflow_index;
CAMELS_BR_data.hfd_mean = hfd_mean;
CAMELS_BR_data.q5 = q5;
CAMELS_BR_data.q95 = q95;
CAMELS_BR_data.high_q_freq = high_q_freq;
CAMELS_BR_data.high_q_dur = high_q_dur;
CAMELS_BR_data.low_q_freq = low_q_freq;
CAMELS_BR_data.low_q_dur = low_q_dur;
CAMELS_BR_data.zero_q_freq = zero_q_freq;

% soil
CAMELS_BR_data.sand_frac = sand_frac;
CAMELS_BR_data.silt_frac = silt_frac;
CAMELS_BR_data.clay_frac = clay_frac;
CAMELS_BR_data.org_carbon_content = org_carbon_content;
CAMELS_BR_data.depth_bedrock = depth_bedrock;
CAMELS_BR_data.water_table_depth = water_table_depth;

% geology
CAMELS_BR_data.geol_1st_class = geol_1st_class;
CAMELS_BR_data.glim_1st_class_frac = glim_1st_class_frac;
CAMELS_BR_data.geol_2nd_class = geol_2nd_class;
CAMELS_BR_data.glim_2nd_class_frac = glim_2nd_class_frac;
CAMELS_BR_data.carbonate_rocks_frac = carbonate_rocks_frac;
CAMELS_BR_data.geol_porosity = geol_porosity;
CAMELS_BR_data.geol_permeability = geol_permeability;

% land cover 
CAMELS_BR_data.crop_frac = crop_frac;
CAMELS_BR_data.crop_mosaic_frac = crop_mosaic_frac;
CAMELS_BR_data.forest_frac = forest_frac;
CAMELS_BR_data.shrub_frac = shrub_frac;
CAMELS_BR_data.grass_frac = grass_frac;
CAMELS_BR_data.bare_frac = bare_frac;
CAMELS_BR_data.imperv_frac = imperv_frac;
CAMELS_BR_data.wet_frac = wet_frac;
CAMELS_BR_data.snow_frac = snow_frac;
CAMELS_BR_data.dom_land_cover = dom_land_cover;
CAMELS_BR_data.dom_land_cover_frac = dom_land_cover_frac;

% human intervention
CAMELS_BR_data.consumptive_use_mm = consumptive_use_mm;
CAMELS_BR_data.consumptive_use = consumptive_use;
CAMELS_BR_data.vol_reservoirs = vol_reservoirs;
CAMELS_BR_data.degree_of_regulation = degree_of_regulation;

% quality check
CAMELS_BR_data.q_quality_control_frac = q_quality_control_frac;
CAMELS_BR_data.q_stream_stage_frac = q_stream_stage_frac;

% hydro-meteorological time series
CAMELS_BR_data.flow_perc_complete = flow_perc_complete;
CAMELS_BR_data.P = P;
CAMELS_BR_data.PET = PET;
CAMELS_BR_data.Q = Q;
CAMELS_BR_data.T = T;

% save file to table
save('./CAMELS_Matlab/Data/CAMELS_BR_data.mat','CAMELS_BR_data')
