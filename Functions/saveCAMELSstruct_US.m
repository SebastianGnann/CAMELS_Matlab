function [CAMELS_US_data] = saveCAMELSstruct_US(save_struct)
%saveCAMELSstruct_US Creates struct file with CAMELS US data.
%   - Loads hydro-meteorological time series and catchment attributes
%   - Timeseries are loaded for the period in which all data are available 
%   - Uses local paths (which contain large CAMELS files)
%   - Adjusts PET data to use standard Priestley-Taylor coefficient of 1.26
%   - Data can be found at: https://ral.ucar.edu/solutions/products/camels
%
%   INPUT
%   save_struct: whether to save struct file or not
%
%   OUTPUT
%   CAMELS_US_data: struct file with CAMELS US data
%
%   References
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
%   Copyright (C) 2021
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 1
    save_struct = false;
end

%% Specify paths
% We have to be above the CAMELS_Matlab directory and CAMELS data should be 
% stored in a folder named CAMELS_US. The following folders are required:
% CAMELS_US/camels_attributes_v2.0/camels_*.txt 
% (7 files; contain catchment attributes)
% CAMELS_US/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/*/*_model_output.txt 
% (18 folders with >1000 files; contain forcing time series)
% CAMELS_US/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2/usgs_streamflow/*/*_streamflow_qc.txt  
% (18 folders with 671 files; contain streamflow time series)

path_catchment_attributes = "CAMELS_US/camels_attributes_v2.0/camels_attributes_v2.0/";
path_modelled_time_series = "CAMELS_US/basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet/"; 
path_observed_time_series = "CAMELS_US/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2/usgs_streamflow/"; 

if ~(exist(path_catchment_attributes) == 7)
    error('Cannot find local path. You can download CAMELS US from https://ral.ucar.edu/solutions/products/camels.')
elseif ~(exist(path_modelled_time_series) == 7)
    error('Cannot find local path. You can download CAMELS US from https://ral.ucar.edu/solutions/products/camels.')
elseif ~(exist(path_observed_time_series) == 7)
    error('Cannot find local path. You can download CAMELS US from https://ral.ucar.edu/solutions/products/camels.')
end

%% Load catchment attributes
% We first load the catchment attribute data which are saved in several txt 
% files.

% topography
% gauge_id;gauge_lat;gauge_lon;elev_mean;slope_mean;area_gages2;area_geospa_fabric
% file_ID_topo = fopen(strcat(path_catchment_attributes,'camels_topo.txt'),'r');
% file_topo = fread(file_ID_topo,'*char');
% file_topo = strrep(file_topo','NA','NaN'); % replace NA with NaN
% camels_topo_data = textscan(file_topo,'%f %f %f %f %f %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_topo);
camels_topo_data = readtable(strcat(path_catchment_attributes,'camels_topo.txt'));

% climate
% gauge_id;p_mean;pet_mean;p_seasonality;frac_snow;aridity;high_prec_freq;high_prec_dur;high_prec_timing;low_prec_freq;low_prec_dur;low_prec_timing
% file_ID_clim = fopen(strcat(path_catchment_attributes,'camels_clim.txt'),'r');
% file_clim = fread(file_ID_clim,'*char');
% file_clim = strrep(file_clim','NA','NaN'); % replace NA with NaN
% camels_climate_data = textscan(file_clim,'%f %f %f %f %f %f %f %f %s %f %f %s',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_clim);
camels_climate_data = readtable(strcat(path_catchment_attributes,'camels_clim.txt'));

% hydrology
% gauge_id;q_mean;runoff_ratio;slope_fdc;baseflow_index;stream_elas;q5;q95;high_q_freq;high_q_dur;low_q_freq;low_q_dur;zero_q_freq;hfd_mean
% file_ID_hydro = fopen(strcat(path_catchment_attributes,'camels_hydro.txt'),'r');
% file_hydro = fread(file_ID_hydro,'*char');
% file_hydro = strrep(file_hydro','NA','NaN'); % replace NA with NaN
% camels_hydro_data = textscan(file_hydro,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_hydro);
camels_hydro_data = readtable(strcat(path_catchment_attributes,'camels_hydro.txt'));

% soils
% gauge_id;soil_depth_pelletier;soil_depth_statsgo;soil_porosity;soil_conductivity;max_water_content;sand_frac;silt_frac;clay_frac;water_frac;organic_frac;other_frac
% file_ID_soil = fopen(strcat(path_catchment_attributes,'camels_soil.txt'),'r');
% file_soil = fread(file_ID_soil,'*char');
% file_soil = strrep(file_soil','NA','NaN'); % replace NA with NaN
% camels_soil_data = textscan(file_soil,'%f %f %f %f %f %f %f %f %f %f %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_soil);
camels_soil_data = readtable(strcat(path_catchment_attributes,'camels_soil.txt'));

% geology
% gauge_id;geol_1st_class;glim_1st_class_frac;geol_2nd_class;glim_2nd_class_frac;carbonate_rocks_frac;geol_porostiy;geol_permeability
% file_ID_geol = fopen(strcat(path_catchment_attributes,'camels_geol.txt'),'r');
% file_geol = fread(file_ID_geol,'*char');
% file_geol = strrep(file_geol','NA','NaN'); % replace NA with NaN
% camels_geol_data = textscan(file_geol,'%f %s %f %s %f %f %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_geol);
camels_geol_data = readtable(strcat(path_catchment_attributes,'camels_geol.txt'));

% vegetation
% gauge_id;frac_forest;lai_max;lai_diff;gvf_max;gvf_diff;dom_land_cover_frac;dom_land_cover;root_depth_50;root_depth_99
% file_ID_vege = fopen(strcat(path_catchment_attributes,'camels_vege.txt'),'r');
% file_vege = fread(file_ID_vege,'*char');
% file_vege = strrep(file_vege','NA','NaN'); % replace NA with NaN
% camels_vege_data = textscan(file_vege,'%f %f %f %f %f %f %f %s %f %f',...
%     'Delimiter',';','headerlines',1);
% fclose(file_ID_vege);
camels_vege_data = readtable(strcat(path_catchment_attributes,'camels_vege.txt'));

% We now add the catchment attributes and metadata to the struct file.

% topography and metadata
CAMELS_US_data.gauge_id = camels_topo_data.gauge_id;
CAMELS_US_data.gauge_lat = camels_topo_data.gauge_lat;
CAMELS_US_data.gauge_lon = camels_topo_data.gauge_lon;
CAMELS_US_data.elev_mean = camels_topo_data.elev_mean;
CAMELS_US_data.slope_mean = camels_topo_data.slope_mean;
CAMELS_US_data.area_gages2 = camels_topo_data.area_gages2;
CAMELS_US_data.area_geospa_fabric = camels_topo_data.area_geospa_fabric;

% climate
CAMELS_US_data.p_mean = camels_climate_data.p_mean;
CAMELS_US_data.pet_mean = camels_climate_data.pet_mean;
CAMELS_US_data.p_seasonality = camels_climate_data.p_seasonality;
CAMELS_US_data.frac_snow = camels_climate_data.frac_snow;
CAMELS_US_data.aridity = camels_climate_data.aridity;
CAMELS_US_data.high_prec_freq = camels_climate_data.high_prec_freq;
CAMELS_US_data.high_prec_dur = camels_climate_data.high_prec_dur;
CAMELS_US_data.high_prec_timing = camels_climate_data.high_prec_timing;
CAMELS_US_data.low_prec_freq = camels_climate_data.low_prec_freq;
CAMELS_US_data.low_prec_dur = camels_climate_data.low_prec_dur;
CAMELS_US_data.low_prec_timing = camels_climate_data.low_prec_timing;

% hydrology
CAMELS_US_data.q_mean = camels_hydro_data.q_mean;
CAMELS_US_data.runoff_ratio = camels_hydro_data.runoff_ratio;
CAMELS_US_data.slope_fdc = camels_hydro_data.slope_fdc;
CAMELS_US_data.baseflow_index = camels_hydro_data.baseflow_index;
CAMELS_US_data.stream_elas = camels_hydro_data.stream_elas;
CAMELS_US_data.q5 = camels_hydro_data.q5;
CAMELS_US_data.q95 = camels_hydro_data.q95;
CAMELS_US_data.high_q_freq = camels_hydro_data.high_q_freq;
CAMELS_US_data.high_q_dur = camels_hydro_data.high_q_dur;
CAMELS_US_data.low_q_freq = camels_hydro_data.low_q_freq;
CAMELS_US_data.low_q_dur = camels_hydro_data.low_q_dur;
CAMELS_US_data.zero_q_freq = camels_hydro_data.zero_q_freq;
CAMELS_US_data.hfd_mean = camels_hydro_data.hfd_mean;

% soil
CAMELS_US_data.soil_depth_pelletier = camels_soil_data.soil_depth_pelletier;
CAMELS_US_data.soil_depth_statsgo = camels_soil_data.soil_depth_statsgo;
CAMELS_US_data.soil_porosity = camels_soil_data.soil_porosity;
CAMELS_US_data.soil_conductivity = camels_soil_data.soil_conductivity;
CAMELS_US_data.max_water_content = camels_soil_data.max_water_content;
CAMELS_US_data.sand_frac = camels_soil_data.sand_frac;
CAMELS_US_data.silt_frac = camels_soil_data.silt_frac;
CAMELS_US_data.clay_frac = camels_soil_data.clay_frac;
CAMELS_US_data.water_frac = camels_soil_data.water_frac;
CAMELS_US_data.organic_frac = camels_soil_data.organic_frac;
CAMELS_US_data.other_frac = camels_soil_data.other_frac;

% geology
CAMELS_US_data.geol_1st_class = camels_geol_data.geol_1st_class;
CAMELS_US_data.glim_1st_class_frac = camels_geol_data.glim_1st_class_frac;
CAMELS_US_data.geol_2nd_class = camels_geol_data.geol_2nd_class;
CAMELS_US_data.glim_2nd_class_frac = camels_geol_data.glim_2nd_class_frac;
CAMELS_US_data.carbonate_rocks_frac = camels_geol_data.carbonate_rocks_frac;
CAMELS_US_data.geol_porosity = camels_geol_data.geol_porostiy;
CAMELS_US_data.geol_permeability = camels_geol_data.geol_permeability;

% vegetation
CAMELS_US_data.frac_forest = camels_vege_data.frac_forest;
CAMELS_US_data.lai_max = camels_vege_data.lai_max;
CAMELS_US_data.lai_diff = camels_vege_data.lai_diff;
CAMELS_US_data.gvf_max = camels_vege_data.gvf_max;
CAMELS_US_data.gvf_diff = camels_vege_data.gvf_diff;
CAMELS_US_data.dom_land_cover_frac = camels_vege_data.dom_land_cover_frac;
CAMELS_US_data.dom_land_cover = camels_vege_data.dom_land_cover;
CAMELS_US_data.root_depth_50 = camels_vege_data.root_depth_50;
CAMELS_US_data.root_depth_99 = camels_vege_data.root_depth_99;

%% Load hydro-meteorological time series
% To extract the time series, we loop over all catchments. We also
% calculate the completeness of the flow records.
flow_perc_complete = NaN(length(CAMELS_US_data.gauge_id),1); % completeness of flow record
P = cell(length(CAMELS_US_data.gauge_id),1); % precipitation
PET = cell(length(CAMELS_US_data.gauge_id),1); % potential evapotranspiration
Q = cell(length(CAMELS_US_data.gauge_id),1); % streamflow
T = cell(length(CAMELS_US_data.gauge_id),1); % temperature

fprintf('Loading catchment data (US)...\n')
for i = 1:length(CAMELS_US_data.gauge_id) % loop over all catchments
    
    if mod(i,100) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,length(CAMELS_US_data.gauge_id))
    end
    
    [P{i}, PET{i}, Q{i}, T{i}] = ...
        loadCatchmentCAMELS_US(CAMELS_US_data.gauge_id(i),...
        path_modelled_time_series,path_observed_time_series,CAMELS_US_data.area_gages2(i));
    flow_perc_complete(i) = 100*(1-sum(isnan(Q{i}(:,2)))./length(Q{i}(:,2)));
    
end

% add hydro-meteorological time series to struct file
CAMELS_US_data.flow_perc_complete = flow_perc_complete;
CAMELS_US_data.P = P;
CAMELS_US_data.PET = PET;
CAMELS_US_data.Q = Q;
CAMELS_US_data.T = T;

% save the struct file
if save_struct
    save('CAMELS_Matlab/Data/CAMELS_US_data.mat','CAMELS_US_data')
end

end