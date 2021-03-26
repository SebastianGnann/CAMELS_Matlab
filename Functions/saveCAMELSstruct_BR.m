function [CAMELS_BR_data] = saveCAMELSstruct_BR(save_struct)
%saveCAMELSstruct_BR Creates struct file with CAMELS BR data.
%   - Loads hydro-meteorological time series and catchment attributes
%   - Timeseries are loaded for the period in which all data are available 
%   - Uses local paths (which contain large CAMELS files)
%   - Loads CHIRPS P and GLEAM PET data
%   - Loads data for 897 selected catchments
%   - Data can be found at: https://zenodo.org/record/3964745
%
%   INPUT
%   save_struct: whether to save struct file or not
%
%   OUTPUT
%   CAMELS_BR_data: struct file with CAMELS BR data
%
%   References
%   Chagas, V.B., Chaffe, P.L., Addor, N., Fan, F.M., Fleischmann, A.S., 
%   Paiva, R.C. and Siqueira, V.A., 2020. CAMELS-BR: hydrometeorological 
%   time series and landscape attributes for 897 catchments in Brazil. 
%   Earth System Science Data, 12(3), pp.2075-2096.
%
%   Copyright (C) 2021
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 1
    save_struct = false;
end


%% Specify paths
% We have to be above the CAMELS_Matlab directory and CAMELS data should be 
% stored in a folder named CAMELS_BR. The following folders are required 
% (each folder contains multiple files):
% CAMELS_BR/1_CAMELS_BR_attributes
% CAMELS_BR/3_CAMELS_BR_streamflow_mm
% CAMELS_BR/4_CAMELS_BR_precipitation_chirps
% CAMELS_BR/9_CAMELS_BR_potential_evapotransp_gleam
% CAMELS_BR/11_CAMELS_BR_temperature_mean_cpc

path_catchment_attributes = "CAMELS_BR/01_CAMELS_BR_attributes/";
path_time_series = "CAMELS_BR/"; 

if ~(exist(path_catchment_attributes) == 7)
    error('Cannot find local path. You can download CAMELS BR from https://zenodo.org/record/3964745.')
elseif ~(exist(path_time_series) == 7)
    error('Cannot find local path. You can download CAMELS BR from https://zenodo.org/record/3964745.')
end

%% Load catchment attributes
% We first load the catchment attribute data which are saved in several txt 
% files.

% location
% gauge_id gauge_name gauge_region gauge_lat gauge_lon
% file_ID_loca = fopen(strcat(path_catchment_attributes,'camels_br_location.txt'),'r');
% file_loca = fread(file_ID_loca,'*char');
% camels_BR_loca_data = textscan(file_loca,'%f %s %s %f %f',...
%     'Delimiter',' ','headerlines',1);
% fclose(file_ID_loca);
camels_BR_loca_data = readtable(strcat(path_catchment_attributes,'camels_br_location.txt'));

% topography
% gauge_id elev_gauge elev_mean slope_mean area
% file_ID_topo = fopen(strcat(path_catchment_attributes,'camels_br_topography.txt'),'r');
% file_topo = fread(file_ID_topo,'*char');
% camels_BR_topo_data = textscan(file_topo,'%f %f %f %f %f',...
%     'Delimiter',' ','headerlines',1);
% fclose(file_ID_topo);
camels_BR_topo_data = readtable(strcat(path_catchment_attributes,'camels_br_topography.txt'));

% climate
% gauge_id p_mean pet_mean aridity p_seasonality frac_snow high_prec_freq high_prec_dur high_prec_timing low_prec_freq low_prec_dur low_prec_timing
% file_ID_im = fread(file_ID_clim,'*char');
% camels_BR_climate_data = textscan(file_clim,'%f %f %f %f %f %f %f %f %s %f %f %s',...
%     'Delimiter',' ','headerlines',1);
% fclose(clim = fopen(strcat(path_catchment_attributes,'camels_br_climate.txt'),'r');
% file_clfile_ID_clim);
camels_BR_climate_data = readtable(strcat(path_catchment_attributes,'camels_br_climate.txt'));

% hydrology
% gauge_id q_mean runoff_ratio stream_elas slope_fdc baseflow_index hfd_mean Q5 Q95 high_q_freq high_q_dur low_q_freq low_q_dur zero_q_freq
% file_ID_hydro = fopen(strcat(path_catchment_attributes,'camels_br_hydrology.txt'),'r');
% file_hydro = fread(file_ID_hydro,'*char');
% camels_BR_hydro_data = textscan(file_hydro,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
%     'Delimiter',' ','headerlines',1);
% fclose(file_ID_hydro);
camels_BR_hydro_data = readtable(strcat(path_catchment_attributes,'camels_br_hydrology.txt'));

% soils
% gauge_id sand_frac silt_frac clay_frac org_carbon_content depth_bedrock water_table_depth
% file_ID_soil = fopen(strcat(path_catchment_attributes,'camels_br_soil.txt'),'r');
% file_soil = fread(file_ID_soil,'*char');
% camels_BR_soil_data = textscan(file_soil,'%f %f %f %f %f %f %f',...
%     'Delimiter',' ','headerlines',1);
% fclose(file_ID_soil);
camels_BR_soil_data = readtable(strcat(path_catchment_attributes,'camels_br_soil.txt'));

% geology
% gauge_id geol_class_1st geol_class_1st_frac geol_class_2nd geol_class_2nd_frac carb_rocks_frac geol_porosity geol_permeability
% file_ID_geol = fopen(strcat(path_catchment_attributes,'camels_br_geology.txt'),'r');
% file_geol = fread(file_ID_geol,'*char');
% camels_BR_geol_data = textscan(file_geol,'%f %s %f %s %f %f %f %f',...
%     'Delimiter',' ','headerlines',1);
% fclose(file_ID_geol);
camels_BR_geol_data = readtable(strcat(path_catchment_attributes,'camels_br_geology.txt'));

% land cover
% gauge_id crop_frac crop_mosaic_frac forest_frac shrub_frac grass_frac bare_frac imperv_frac wet_frac snow_frac dom_land_cover dom_land_cover_frac
% file_ID_land = fopen(strcat(path_catchment_attributes,'camels_br_land_cover.txt'),'r');
% file_land = fread(file_ID_land,'*char');
% camels_BR_land_data = textscan(file_land,'%f %f %f %f %f %f %f %f %f %f %s %f',...
%     'Delimiter',' ','headerlines',1);
% fclose(file_ID_land);
camels_BR_land_data = readtable(strcat(path_catchment_attributes,'camels_br_land_cover.txt'));

% human intervention
% gauge_id consumptive_use_mm consumptive_use vol_reservoirs degree_of_regulation
% file_ID_huma = fopen(strcat(path_catchment_attributes,'camels_br_human_intervention.txt'),'r');
% file_huma = fread(file_ID_huma,'*char');
% camels_BR_huma_data = textscan(file_huma,'%f %f %f %f %f',...
%     'Delimiter',' ','headerlines',1);
% fclose(file_ID_huma);
camels_BR_huma_data = readtable(strcat(path_catchment_attributes,'camels_br_human_intervention.txt'));

% quality check 
% gauge_id q_quality_control_frac q_stream_stage_frac
% file_ID_qual = fopen(strcat(path_catchment_attributes,'camels_br_quality_check.txt'),'r');
% file_qual = fread(file_ID_qual,'*char');
% camels_BR_qual_data = textscan(file_qual,'%f %f %f',...
%     'Delimiter',' ','headerlines',1);
% fclose(file_ID_qual);
camels_BR_qual_data = readtable(strcat(path_catchment_attributes,'camels_br_quality_check.txt'));

% We now add the catchment attributes and metadata to the struct file.

CAMELS_BR_data.gauge_id = camels_BR_climate_data.gauge_id;

% location
gauge_id_tmp = camels_BR_loca_data.gauge_id;
gauge_name_tmp = camels_BR_loca_data.gauge_name;
gauge_region_tmp = camels_BR_loca_data.gauge_region;
gauge_lat_tmp = camels_BR_loca_data.gauge_lat;
gauge_lon_tmp = camels_BR_loca_data.gauge_lon;
% location is provided for more than the 897 selected catchments, so we
% need to pick only the location of the 897 catchments of interest
CAMELS_BR_data.gauge_name = strings(size(CAMELS_BR_data.gauge_id));
CAMELS_BR_data.gauge_region = strings(size(CAMELS_BR_data.gauge_id));
CAMELS_BR_data.gauge_lat = NaN(size(CAMELS_BR_data.gauge_id));
CAMELS_BR_data.gauge_lon  = NaN(size(CAMELS_BR_data.gauge_id));
for i = 1:length(CAMELS_BR_data.gauge_name)
    index = find(gauge_id_tmp==CAMELS_BR_data.gauge_id(i));
    CAMELS_BR_data.gauge_name(i) = gauge_name_tmp(index);
    CAMELS_BR_data.gauge_region(i) = gauge_region_tmp(index);
    CAMELS_BR_data.gauge_lat(i) = gauge_lat_tmp(index);
    CAMELS_BR_data.gauge_lon(i) = gauge_lon_tmp(index);
end

% topography
CAMELS_BR_data.elev_gauge = camels_BR_topo_data.elev_gauge;
CAMELS_BR_data.elev_mean = camels_BR_topo_data.elev_mean;
CAMELS_BR_data.slope_mean = camels_BR_topo_data.slope_mean;
CAMELS_BR_data.area = camels_BR_topo_data.area;

% climate
CAMELS_BR_data.p_mean = camels_BR_climate_data.p_mean;
CAMELS_BR_data.pet_mean = camels_BR_climate_data.pet_mean;
CAMELS_BR_data.aridity = camels_BR_climate_data.aridity;
CAMELS_BR_data.p_seasonality = camels_BR_climate_data.p_seasonality;
CAMELS_BR_data.frac_snow = camels_BR_climate_data.frac_snow;
CAMELS_BR_data.high_prec_freq = camels_BR_climate_data.high_prec_freq;
CAMELS_BR_data.high_prec_dur = camels_BR_climate_data.high_prec_dur;
CAMELS_BR_data.high_prec_timing = camels_BR_climate_data.high_prec_timing;
CAMELS_BR_data.low_prec_freq = camels_BR_climate_data.low_prec_freq;
CAMELS_BR_data.low_prec_dur = camels_BR_climate_data.low_prec_dur;
CAMELS_BR_data.low_prec_timing = camels_BR_climate_data.low_prec_timing;

% hydrology 
CAMELS_BR_data.q_mean = camels_BR_hydro_data.q_mean;
CAMELS_BR_data.runoff_ratio = camels_BR_hydro_data.runoff_ratio;
CAMELS_BR_data.stream_elas = camels_BR_hydro_data.stream_elas;
CAMELS_BR_data.slope_fdc = camels_BR_hydro_data.slope_fdc;
CAMELS_BR_data.baseflow_index = camels_BR_hydro_data.baseflow_index;
CAMELS_BR_data.hfd_mean = camels_BR_hydro_data.hfd_mean;
CAMELS_BR_data.q5 = camels_BR_hydro_data.Q5;
CAMELS_BR_data.q95 = camels_BR_hydro_data.Q95;
CAMELS_BR_data.high_q_freq = camels_BR_hydro_data.high_q_freq;
CAMELS_BR_data.high_q_dur = camels_BR_hydro_data.high_q_freq;
CAMELS_BR_data.low_q_freq = camels_BR_hydro_data.low_q_freq;
CAMELS_BR_data.low_q_dur = camels_BR_hydro_data.low_q_dur;
CAMELS_BR_data.zero_q_freq = camels_BR_hydro_data.zero_q_freq;

% soil
CAMELS_BR_data.sand_perc = camels_BR_soil_data.sand_perc;
CAMELS_BR_data.silt_perc = camels_BR_soil_data.silt_perc;
CAMELS_BR_data.clay_perc = camels_BR_soil_data.clay_perc;
CAMELS_BR_data.org_carbon_content = camels_BR_soil_data.org_carbon_content;
CAMELS_BR_data.bedrock_depth = camels_BR_soil_data.bedrock_depth;
CAMELS_BR_data.water_table_depth = camels_BR_soil_data.water_table_depth;

% geology
CAMELS_BR_data.geol_class_1st = camels_BR_geol_data.geol_class_1st;
CAMELS_BR_data.geol_class_1st_perc = camels_BR_geol_data.geol_class_1st_perc;
CAMELS_BR_data.geol_class_2nd = camels_BR_geol_data.geol_class_2nd;
CAMELS_BR_data.geol_class_2nd_perc = camels_BR_geol_data.geol_class_2nd_perc;
CAMELS_BR_data.carb_rocks_perc = camels_BR_geol_data.carb_rocks_perc;
CAMELS_BR_data.geol_porosity = camels_BR_geol_data.geol_porosity;
CAMELS_BR_data.geol_permeability = camels_BR_geol_data.geol_permeability;

% land cover 
CAMELS_BR_data.crop_perc = camels_BR_land_data.crop_perc;
CAMELS_BR_data.crop_mosaic_perc = camels_BR_land_data.crop_mosaic_perc;
CAMELS_BR_data.forest_perc = camels_BR_land_data.forest_perc;
CAMELS_BR_data.shrub_perc = camels_BR_land_data.shrub_perc;
CAMELS_BR_data.grass_perc = camels_BR_land_data.grass_perc;
CAMELS_BR_data.barren_perc = camels_BR_land_data.barren_perc;
CAMELS_BR_data.imperv_perc = camels_BR_land_data.imperv_perc;
CAMELS_BR_data.wet_perc = camels_BR_land_data.wet_perc;
CAMELS_BR_data.snow_perc = camels_BR_land_data.snow_perc;
CAMELS_BR_data.dom_land_cover = camels_BR_land_data.dom_land_cover;
CAMELS_BR_data.dom_land_cover_perc = camels_BR_land_data.dom_land_cover_perc;

% human intervention
CAMELS_BR_data.consumptive_use = camels_BR_huma_data.consumptive_use;
CAMELS_BR_data.consumptive_use_perc = camels_BR_huma_data.consumptive_use_perc;
CAMELS_BR_data.reservoirs_vol = camels_BR_huma_data.reservoirs_vol;
CAMELS_BR_data.regulation_degree = camels_BR_huma_data.regulation_degree;

% quality check
CAMELS_BR_data.q_quality_control_frac = camels_BR_qual_data.q_quality_control_perc;
CAMELS_BR_data.q_stream_stage_frac = camels_BR_qual_data.q_stream_stage_perc;

%% Load hydro-meteorological time series
% To extract the time series, we loop over all catchments. We also
% calculate the completeness of the flow records.
flow_perc_complete = NaN(length(CAMELS_BR_data.gauge_id),1);
P = cell(length(CAMELS_BR_data.gauge_id),1); % precipitation
PET = cell(length(CAMELS_BR_data.gauge_id),1); % potential evapotranspiration
Q = cell(length(CAMELS_BR_data.gauge_id),1); % streamflow
T = cell(length(CAMELS_BR_data.gauge_id),1); % temperature

fprintf('Loading catchment data (BR)...\n')
for i = 1:length(CAMELS_BR_data.gauge_id) % loop over all catchments
    
    if mod(i,100) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,length(CAMELS_BR_data.gauge_id))
    end
    
    [P{i}, PET{i}, Q{i}, T{i}] = ...
        loadCatchmentCAMELS_BR(CAMELS_BR_data.gauge_id(i),path_time_series);
    flow_perc_complete(i) = 100*(1-sum(isnan(Q{i}(:,2)))./length(Q{i}(:,2)));
    
end

% add hydro-meteorological time series to struct file
CAMELS_BR_data.flow_perc_complete = flow_perc_complete;
CAMELS_BR_data.P = P;
CAMELS_BR_data.PET = PET;
CAMELS_BR_data.Q = Q;
CAMELS_BR_data.T = T;

% save the struct file
if save_struct
    save('CAMELS_Matlab/Data/CAMELS_BR_data.mat','-struct','CAMELS_BR_data')
end

end
