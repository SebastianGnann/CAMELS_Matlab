%% Data preparation for CAMELS-GB data:
%   - Create structure with CAMELS GB catchment data and attributes
%   - Uses local paths (large CAMELS-GB files)
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

clc
clear all
close all

%% add directories for functions to path

%% load CAMELS data

% Before you can run the code you need to extract CAMELS data and store
% them in the correct folder. The following folders are required:
%   - camels_attributes_v2.0
%   - basin_timeseries_v1p2_modelOutput_daymet (used for time series)

disp('This function uses local paths. Change to local path where CAMELS-GB data are stored.')

path_camelsGB = "C:\Users\sg16200\Local Documents\CAMELS_GB\"; % change to your local path
path_time_series = strcat(path_camelsGB,'CAMELS_GB_Timeseries\'); % change to your local path
path_catchment_attributes = strcat(path_camelsGB,'CAMELSGB_Catchment_Attributes.xlsx'); % change to your local path

if exist(path_camelsGB) == 7
    addpath(genpath(path_camelsGB)); % not actually needed since the files are just loaded
else
    error('Cannot find local path. You can download CAMELS-GB from <tba>.')
end

%% load catchment attributes
% load data from files

% topography
% gauge_id	gauge_name	gauge_lat	gauge_lon	gauge_easting	gauge_northing	gauge_elev	area	dpsbar	elev_mean	elev_min	elev_10	elev_50	elev_90	elev_max
[camelsGB_topo_data,camelsGB_topo_data_str] = xlsread(path_catchment_attributes,1);

% climatic indices
% gauge_id	p_mean	pet_mean	aridity	p_seasonality	frac_snow	high_prec_freq	high_prec_dur	high_prec_timing	low_prec_freq	low_prec_dur	low_prec_timing
[camelsGB_climate_data,camelsGB_climate_data_str] = xlsread(path_catchment_attributes,2);

% hydrology
% gauge_id	q_mean	runoff_ratio	stream_elas	slope_fdc	baseflow_index	baseflow_index_ceh	hfd_mean	Q5	Q95	high_q_freq	high_q_dur	low_q_freq	low_q_dur	zero_q_freq
[camelsGB_hydro_data,camelsGB_hydro_data_str] = xlsread(path_catchment_attributes,3);

% land cover
% gauge_id	dwood_perc	ewood_perc	grass_perc	shrub_perc	crop_perc	urban_perc	inwater_perc	bares_perc	dom_land_cover
[camelsGB_land_data,camelsGB_land_data_str] = xlsread(path_catchment_attributes,4);

% soils
% gauge_id	sand_perc	sand_perc_missing	silt_perc	silt_perc_missing	clay_perc	clay_perc_missing	organic_perc	organic_perc_missing	bulkdens	bulkdens_missing	tawc	tawc_missing	porosity_cosby	porosity_cosby_missing	porosity_hypres	porosity_hypres_missing	conductivity_cosby	conductivity_cosby_missing	conductivity_hypres	conductivity_hypres_missing	root_depth	root_depth_missing	soil_depth_pelletier	soil_depth_pelletier_missing
[camelsGB_soil_data,camelsGB_soil_data_str] = xlsread(path_catchment_attributes,5);

% hydrogeology
% gauge_id	inter_high_perc	inter_mod_perc	inter_low_perc	frac_high_perc	frac_mod_perc	frac_low_perc	no_gw_perc	low_nsig_perc	nsig_low_perc
[camelsGB_geology_data,camelsGB_geology_data_str] = xlsread(path_catchment_attributes,6);

% hydrometry data
% gauge_id	station_type	flow_period_start	flow_period_end	flow_perc_complete	bankfull_flow	structurefull_flow	num_rcurves	num_gaugings	max_gauging_flow	max_gauging_level	max_gauging_date	min_gauging_flow	min_gauging_level	min_gauging_date	q1_uncert_upper	q1_uncert_lower	q5_uncert_upper	q5_uncert_lower	q25_uncert_upper	q25_uncert_lower	q50_uncert_upper	q50_uncert_lower	q95_uncert_upper	q95_uncert_lower	q99_uncert_upper	q99_uncert_lower
[camelsGB_hydrometry_data,camelsGB_hydrometry_data_str] = xlsread(path_catchment_attributes,7);

% human influences
% gauge_id	benchmark_catch	num_reservoir	reservoir_cap	reservoir_he	reservoir_nav	reservoir_drain	reservoir_wr	reservoir_fs	reservoir_env	reservoir_nousedata	reservoir_year_first	reservoir_year_last
[camelsGB_human_data,camelsGB_human_data_str] = xlsread(path_catchment_attributes,8);

%% extract data from matrices

% topography
gauge_id = camelsGB_topo_data(:,1);
gauge_name = camelsGB_topo_data_str(2:end,2); %
gauge_lat = camelsGB_topo_data(:,3);
gauge_lon = camelsGB_topo_data(:,4);
gauge_easting = camelsGB_topo_data(:,5);
gauge_northing = camelsGB_topo_data(:,6);
gauge_elev = camelsGB_topo_data(:,7);
area = camelsGB_topo_data(:,8);
dpsbar = camelsGB_topo_data(:,9);
elev_mean = camelsGB_topo_data(:,10);
elev_min = camelsGB_topo_data(:,11);
elev_10 = camelsGB_topo_data(:,12);
elev_50 = camelsGB_topo_data(:,13);
elev_90 = camelsGB_topo_data(:,14);
elev_max = camelsGB_topo_data(:,15);

% climatic indices
p_mean = camelsGB_climate_data(:,2);
pet_mean = camelsGB_climate_data(:,3);
aridity = camelsGB_climate_data(:,4);
p_seasonality = camelsGB_climate_data(:,5);
frac_snow = camelsGB_climate_data(:,6);
high_prec_freq = camelsGB_climate_data(:,7);
high_prec_dur = camelsGB_climate_data(:,8);
high_prec_timing = camelsGB_climate_data_str(2:end,9);
low_prec_freq = camelsGB_climate_data(:,10);
low_prec_dur = camelsGB_climate_data(:,11);
low_prec_timing = camelsGB_climate_data_str(2:end,12);

% hydrological signatures
q_mean = camelsGB_hydro_data(:,2);
runoff_ratio = camelsGB_hydro_data(:,3);
stream_elas = camelsGB_hydro_data(:,4);
slope_fdc = camelsGB_hydro_data(:,5);
baseflow_index = camelsGB_hydro_data(:,6);
baseflow_index_ceh = camelsGB_hydro_data(:,7);
hfd_mean = camelsGB_hydro_data(:,8);
q5 = camelsGB_hydro_data(:,9);
q95 = camelsGB_hydro_data(:,10);
high_q_freq = camelsGB_hydro_data(:,11);
high_q_dur = camelsGB_hydro_data(:,12);
low_q_freq = camelsGB_hydro_data(:,13);
low_q_dur = camelsGB_hydro_data(:,14);
zero_q_freq = camelsGB_hydro_data(:,15);

% land cover
dwood_perc = camelsGB_land_data(:,2);
ewood_perc = camelsGB_land_data(:,3);
grass_perc = camelsGB_land_data(:,4);
shrub_perc = camelsGB_land_data(:,5);
crop_perc = camelsGB_land_data(:,6);
urban_perc = camelsGB_land_data(:,7);
inwater_perc = camelsGB_land_data(:,8);
bares_perc = camelsGB_land_data(:,9);
dom_land_cover = camelsGB_land_data_str(2:end,10);

% soils
sand_perc = camelsGB_soil_data(:,2);
sand_perc_missing = camelsGB_soil_data(:,3);
silt_perc = camelsGB_soil_data(:,4);
silt_perc_missing = camelsGB_soil_data(:,5);
clay_perc = camelsGB_soil_data(:,6);
clay_perc_missing = camelsGB_soil_data(:,7);
organic_perc = camelsGB_soil_data(:,8);
organic_perc_missing = camelsGB_soil_data(:,9);
bulkdens = camelsGB_soil_data(:,10);
bulkdens_missing = camelsGB_soil_data(:,11);
tawc = camelsGB_soil_data(:,12);
tawc_missing = camelsGB_soil_data(:,13);
porosity_cosby = camelsGB_soil_data(:,14);
porosity_cosby_missing = camelsGB_soil_data(:,15);
porosity_hypres = camelsGB_soil_data(:,16);
porosity_hypres_missing = camelsGB_soil_data(:,17);
conductivity_cosby = camelsGB_soil_data(:,18);
conductivity_cosby_missing = camelsGB_soil_data(:,19);
conductivity_hypres = camelsGB_soil_data(:,20);
conductivity_hypres_missing = camelsGB_soil_data(:,21);
root_depth = camelsGB_soil_data(:,22);
root_depth_missing = camelsGB_soil_data(:,23);
soil_depth_pelletier = camelsGB_soil_data(:,24);
soil_depth_pelletier_missing = camelsGB_soil_data(:,25);

% hydrogeology
inter_high_perc = camelsGB_geology_data(:,2);
inter_mod_perc = camelsGB_geology_data(:,3);
inter_low_perc = camelsGB_geology_data(:,4);
frac_high_perc = camelsGB_geology_data(:,5);
frac_mod_perc = camelsGB_geology_data(:,6);
frac_low_perc = camelsGB_geology_data(:,7);
no_gw_perc = camelsGB_geology_data(:,8);
low_nsig_perc = camelsGB_geology_data(:,9);
nsig_low_perc = camelsGB_geology_data(:,10);

% hydrometry
station_type = camelsGB_hydrometry_data_str(2:end,2);
flow_period_start = camelsGB_hydrometry_data_str(2:end,3);
flow_period_end = camelsGB_hydrometry_data_str(2:end,4);
flow_perc_complete = camelsGB_hydrometry_data(:,5);
bankfull_flow = camelsGB_hydrometry_data(:,6);
structurefull_flow = camelsGB_hydrometry_data(:,7);
num_rcurves = camelsGB_hydrometry_data(:,8);
num_gaugings = camelsGB_hydrometry_data(:,9);
max_gauging_flow = camelsGB_hydrometry_data(:,10);
max_gauging_level = camelsGB_hydrometry_data(:,11);
max_gauging_date = camelsGB_hydrometry_data_str(2:end,12);
min_gauging_flow = camelsGB_hydrometry_data(:,13);
min_gauging_level = camelsGB_hydrometry_data(:,14);
min_gauging_date = camelsGB_hydrometry_data_str(2:end,15);
q1_uncert_upper = camelsGB_hydrometry_data(:,16);
q1_uncert_lower = camelsGB_hydrometry_data(:,17);
q5_uncert_upper = camelsGB_hydrometry_data(:,18);
q5_uncert_lower = camelsGB_hydrometry_data(:,19);
q25_uncert_upper = camelsGB_hydrometry_data(:,20);
q25_uncert_lower = camelsGB_hydrometry_data(:,21);
q50_uncert_upper = camelsGB_hydrometry_data(:,22);
q50_uncert_lower = camelsGB_hydrometry_data(:,23);
q95_uncert_upper = camelsGB_hydrometry_data(:,24);
q95_uncert_lower = camelsGB_hydrometry_data(:,25);
q99_uncert_upper = camelsGB_hydrometry_data(:,26);
q99_uncert_lower = camelsGB_hydrometry_data(:,27);

% human influences
benchmark_catch = camelsGB_human_data_str(2:end,2);
isBenchmark = zeros(size(benchmark_catch));
isBenchmark(strcmp('Y',benchmark_catch)) = 1;
isBenchmark = boolean(isBenchmark);
num_reservoir = camelsGB_human_data(:,3);
reservoir_cap = camelsGB_human_data(:,4);
reservoir_he = camelsGB_human_data(:,5);
reservoir_nav = camelsGB_human_data(:,6);
reservoir_drain = camelsGB_human_data(:,7);
reservoir_wr = camelsGB_human_data(:,8);
reservoir_fs = camelsGB_human_data(:,9);
reservoir_env = camelsGB_human_data(:,10);
reservoir_nousedata = camelsGB_human_data(:,11);
reservoir_year_first = camelsGB_human_data(:,12);
reservoir_year_last = camelsGB_human_data(:,13);

%% load time series
% intialise arrays to store some metrics and hydrological signatures
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
    [P{i}, PET{i}, Q{i}, T{i}] = loadCatchment_ID_CAMELS_GB(ID,path_time_series);
    loadSuccesful = true;
    
    if sum(isnan(P(:,2))) > 0
        disp('NaN P');
    elseif sum(isnan(PET(:,2))) > 0
        disp('NaN PET');
    elseif sum(isnan(Q(:,2))) > 0
        disp('NaN Q');
    end
    
end

%% create structure

% topography
CAMELS_GB_data.gauge_id = gauge_id;
CAMELS_GB_data.gauge_name = gauge_name;
CAMELS_GB_data.gauge_lat = gauge_lat;
CAMELS_GB_data.gauge_lon = gauge_lon;
CAMELS_GB_data.gauge_easting = gauge_easting;
CAMELS_GB_data.gauge_northing = gauge_northing;
CAMELS_GB_data.gauge_elev = gauge_elev;
CAMELS_GB_data.area = area;
CAMELS_GB_data.dpsbar = dpsbar;
CAMELS_GB_data.elev_mean = elev_mean;
CAMELS_GB_data.elev_min = elev_min;
CAMELS_GB_data.elev_10 = elev_10;
CAMELS_GB_data.elev_50 = elev_50;
CAMELS_GB_data.elev_90 = elev_90;
CAMELS_GB_data.elev_max = elev_max;

% climatic indices
CAMELS_GB_data.p_mean = p_mean;
CAMELS_GB_data.pet_mean = pet_mean;
CAMELS_GB_data.aridity = aridity;
CAMELS_GB_data.p_seasonality = p_seasonality;
CAMELS_GB_data.frac_snow = frac_snow;
CAMELS_GB_data.high_prec_freq = high_prec_freq;
CAMELS_GB_data.high_prec_dur = high_prec_dur;
CAMELS_GB_data.high_prec_timing = high_prec_timing;
CAMELS_GB_data.low_prec_freq = low_prec_freq;
CAMELS_GB_data.low_prec_dur = low_prec_dur;
CAMELS_GB_data.low_prec_timing = low_prec_timing;

% hydrological signatures
CAMELS_GB_data.q_mean = q_mean;
CAMELS_GB_data.runoff_ratio = runoff_ratio;
CAMELS_GB_data.stream_elas = stream_elas;
CAMELS_GB_data.slope_fdc = slope_fdc;
CAMELS_GB_data.baseflow_index = baseflow_index;
CAMELS_GB_data.baseflow_index_ceh = baseflow_index_ceh;
CAMELS_GB_data.hfd_mean = hfd_mean;
CAMELS_GB_data.q5 = q5;
CAMELS_GB_data.q95 = q95;
CAMELS_GB_data.high_q_freq = high_q_freq;
CAMELS_GB_data.high_q_dur = high_q_dur;
CAMELS_GB_data.low_q_freq = low_q_freq;
CAMELS_GB_data.low_q_dur = low_q_dur;
CAMELS_GB_data.zero_q_freq = zero_q_freq;

% land cover
CAMELS_GB_data.dwood_perc = dwood_perc;
CAMELS_GB_data.ewood_perc = ewood_perc;
CAMELS_GB_data.grass_perc = grass_perc;
CAMELS_GB_data.shrub_perc = shrub_perc;
CAMELS_GB_data.crop_perc = crop_perc;
CAMELS_GB_data.urban_perc = urban_perc;
CAMELS_GB_data.inwater_perc = inwater_perc;
CAMELS_GB_data.bares_perc = bares_perc;
CAMELS_GB_data.dom_land_cover = dom_land_cover;

% soils
CAMELS_GB_data.sand_perc = sand_perc;
CAMELS_GB_data.sand_perc_missing = sand_perc_missing;
CAMELS_GB_data.silt_perc = silt_perc;
CAMELS_GB_data.silt_perc_missing = silt_perc_missing;
CAMELS_GB_data.clay_perc = clay_perc;
CAMELS_GB_data.clay_perc_missing = clay_perc_missing;
CAMELS_GB_data.organic_perc = organic_perc;
CAMELS_GB_data.organic_perc_missing = organic_perc_missing;
CAMELS_GB_data.bulkdens = bulkdens;
CAMELS_GB_data.bulkdens_missing = bulkdens_missing;
CAMELS_GB_data.tawc = tawc;
CAMELS_GB_data.tawc_missing = tawc_missing;
CAMELS_GB_data.porosity_cosby = porosity_cosby;
CAMELS_GB_data.porosity_cosby_missing = porosity_cosby_missing;
CAMELS_GB_data.porosity_hypres = porosity_hypres;
CAMELS_GB_data.porosity_hypres_missing = porosity_hypres_missing;
CAMELS_GB_data.conductivity_cosby = conductivity_cosby;
CAMELS_GB_data.conductivity_cosby_missing = conductivity_cosby_missing;
CAMELS_GB_data.conductivity_hypres = conductivity_hypres;
CAMELS_GB_data.conductivity_hypres_missing = conductivity_hypres_missing;
CAMELS_GB_data.root_depth = root_depth;
CAMELS_GB_data.root_depth_missing = root_depth_missing;
CAMELS_GB_data.soil_depth_pelletier = soil_depth_pelletier;
CAMELS_GB_data.soil_depth_pelletier_missing = soil_depth_pelletier_missing;

% hydrogeology
CAMELS_GB_data.inter_high_perc = inter_high_perc;
CAMELS_GB_data.inter_mod_perc = inter_mod_perc;
CAMELS_GB_data.inter_low_perc = inter_low_perc;
CAMELS_GB_data.frac_high_perc = frac_high_perc;
CAMELS_GB_data.frac_mod_perc = frac_mod_perc;
CAMELS_GB_data.frac_low_perc = frac_low_perc;
CAMELS_GB_data.no_gw_perc = no_gw_perc;
CAMELS_GB_data.low_nsig_perc = low_nsig_perc;
CAMELS_GB_data.low_nsig_perc = low_nsig_perc;

% hydrometry
CAMELS_GB_data.station_type = station_type;
CAMELS_GB_data.flow_period_start = flow_period_start;
CAMELS_GB_data.flow_period_end = flow_period_end;
CAMELS_GB_data.flow_perc_complete = flow_perc_complete;
CAMELS_GB_data.bankfull_flow = bankfull_flow;
CAMELS_GB_data.structurefull_flow = structurefull_flow;
CAMELS_GB_data.num_rcurves = num_rcurves;
CAMELS_GB_data.num_gaugings = num_gaugings;
CAMELS_GB_data.max_gauging_flow = max_gauging_flow;
CAMELS_GB_data.max_gauging_level = max_gauging_level;
CAMELS_GB_data.max_gauging_date = max_gauging_date;
CAMELS_GB_data.min_gauging_flow = min_gauging_flow;
CAMELS_GB_data.min_gauging_level = min_gauging_level;
CAMELS_GB_data.min_gauging_date = min_gauging_date;
CAMELS_GB_data.q1_uncert_upper = q1_uncert_upper;
CAMELS_GB_data.q1_uncert_lower = q1_uncert_lower;
CAMELS_GB_data.q5_uncert_upper = q5_uncert_upper;
CAMELS_GB_data.q5_uncert_lower = q5_uncert_lower;
CAMELS_GB_data.q25_uncert_upper = q25_uncert_upper;
CAMELS_GB_data.q25_uncert_lower = q25_uncert_lower;
CAMELS_GB_data.q50_uncert_upper = q50_uncert_upper;
CAMELS_GB_data.q50_uncert_lower = q50_uncert_lower;
CAMELS_GB_data.q95_uncert_upper = q95_uncert_upper;
CAMELS_GB_data.q95_uncert_lower = q95_uncert_lower;
CAMELS_GB_data.q99_uncert_upper = q99_uncert_upper;
CAMELS_GB_data.q99_uncert_lower = q99_uncert_lower;

% human influences
CAMELS_GB_data.benchmark_catch = benchmark_catch;
CAMELS_GB_data.isBenchmark = isBenchmark;
CAMELS_GB_data.num_reservoir = num_reservoir;
CAMELS_GB_data.reservoir_cap = reservoir_cap;
CAMELS_GB_data.reservoir_he = reservoir_he;
CAMELS_GB_data.reservoir_nav = reservoir_nav;
CAMELS_GB_data.reservoir_drain = reservoir_drain;
CAMELS_GB_data.reservoir_wr = reservoir_wr;
CAMELS_GB_data.reservoir_fs = reservoir_fs;
CAMELS_GB_data.reservoir_env = reservoir_env;
CAMELS_GB_data.reservoir_nousedata = reservoir_nousedata;
CAMELS_GB_data.reservoir_year_first = reservoir_year_first;
CAMELS_GB_data.reservoir_year_last = reservoir_year_last;

% hydro-meteorological time series
CAMELS_GB_data.P = P;
CAMELS_GB_data.PET= PET;
CAMELS_GB_data.Q = Q;
CAMELS_GB_data.T = T;

% save file to table
save('./CAMELS_Matlab/Data/CAMELS_GB_data.mat','CAMELS_GB_data')

