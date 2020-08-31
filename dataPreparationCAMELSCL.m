%% Data preparation for CAMELS-CL data:
%   - Create struct file with CAMELS-CL time series and attributes
%   - Uses local paths (large CAMELS-CL files)
%   - Uses MSWEP P and HARGREAVES PET data for time series
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
%   \CAMELS_CL
%       \1_CAMELScl_attributes
%       \3_CAMELScl_streamflow_mm
%       \5_CAMELScl_precip_chirps
%       \6_CAMELScl_precip_mswep
%       \10_CAMELScl_tmean_cr2met
%       \11_CAMELScl_pet_8d_modis
%       \12_CAMELScl_pet_hargreaves

disp('This function uses local paths. Change to local path where CAMELS-CL data are stored.')

path_catchment_attributes = "C:\Users\sg16200\Local Documents\CAMELS_CL\1_CAMELScl_attributes\1_CAMELScl_attributes.txt"; % change to your local path
path_time_series = "C:\Users\sg16200\Local Documents\CAMELS_CL\"; % change to your local path

if ~exist(path_time_series) == 7
    error('Cannot find local path. You can download CAMELS-CL from https://doi.pangaea.de/10.1594/PANGAEA.894885?format=html#download.')
end

%% load catchment attributes
% load data from files

% all attributes
camels_CL_data = readtable(path_catchment_attributes);
% replace NA with NaN
for i = 1:size(camels_CL_data,2)
    char1 = 'NA';
    char2 = 'NaN';
    eval(['camels_CL_data.Var',num2str(i),...
        '(strcmp(camels_CL_data.Var',num2str(i),',char1)) = {char2};'])
end

%% extract data from matrices
% topography and metadata
gauge_id = cellfun(@str2num,table2array(camels_CL_data(1,2:end))');
gauge_name = table2array(camels_CL_data(2,2:end))';
gauge_lat = cellfun(@str2num,table2array(camels_CL_data(3,2:end))');
gauge_lon = cellfun(@str2num,table2array(camels_CL_data(4,2:end))');
% record_period_start
% record_period_end
% n_obs
area = cellfun(@str2num,table2array(camels_CL_data(8,2:end))');
elev_gauge = cellfun(@str2num,table2array(camels_CL_data(9,2:end))');
elev_mean = cellfun(@str2num,table2array(camels_CL_data(10,2:end))');
elev_med = cellfun(@str2num,table2array(camels_CL_data(11,2:end))');
elev_max = cellfun(@str2num,table2array(camels_CL_data(12,2:end))');
elev_min = cellfun(@str2num,table2array(camels_CL_data(13,2:end))');
slope_mean = cellfun(@str2num,table2array(camels_CL_data(14,2:end))');
nested_inner = cellfun(@str2num,table2array(camels_CL_data(15,2:end))');
nested_outer = cellfun(@str2num,table2array(camels_CL_data(16,2:end))');
location_type = table2array(camels_CL_data(17,2:end))';
% geology
geol_class_1st = table2array(camels_CL_data(18,2:end))';
geol_class_1st_frac = cellfun(@str2num,table2array(camels_CL_data(19,2:end))');
geol_class_2nd = table2array(camels_CL_data(1,20:end))';
geol_class_2nd_frac = cellfun(@str2num,table2array(camels_CL_data(21,2:end))');
carb_rocks_frac = cellfun(@str2num,table2array(camels_CL_data(22,2:end))');
% land use
crop_frac = cellfun(@str2num,table2array(camels_CL_data(23,2:end))');
nf_frac = cellfun(@str2num,table2array(camels_CL_data(24,2:end))');
fp_frac = cellfun(@str2num,table2array(camels_CL_data(25,2:end))');
grass_frac = cellfun(@str2num,table2array(camels_CL_data(26,2:end))');
shrub_frac = cellfun(@str2num,table2array(camels_CL_data(27,2:end))');
wet_frac = cellfun(@str2num,table2array(camels_CL_data(28,2:end))');
imp_frac = cellfun(@str2num,table2array(camels_CL_data(29,2:end))');
lc_barren = cellfun(@str2num,table2array(camels_CL_data(30,2:end))');
snow_frac = cellfun(@str2num,table2array(camels_CL_data(31,2:end))');
lc_glacier = cellfun(@str2num,table2array(camels_CL_data(32,2:end))');
fp_nf_index = cellfun(@str2num,table2array(camels_CL_data(33,2:end))');
forest_frac = cellfun(@str2num,table2array(camels_CL_data(34,2:end))');
dom_land_cover = table2array(camels_CL_data(35,2:end))';
dom_land_cover_frac = cellfun(@str2num,table2array(camels_CL_data(36,2:end))');
land_cover_missing = cellfun(@str2num,table2array(camels_CL_data(37,2:end))');
% climate
p_mean_cr2met = cellfun(@str2num,table2array(camels_CL_data(38,2:end))');
p_mean_chirps = cellfun(@str2num,table2array(camels_CL_data(39,2:end))');
p_mean_mswep = cellfun(@str2num,table2array(camels_CL_data(40,2:end))');
p_mean_tmpa = cellfun(@str2num,table2array(camels_CL_data(41,2:end))');
pet_mean = cellfun(@str2num,table2array(camels_CL_data(42,2:end))');
aridity_cr2met = cellfun(@str2num,table2array(camels_CL_data(43,2:end))');
aridity_chirps = cellfun(@str2num,table2array(camels_CL_data(44,2:end))');
aridity_mswep = cellfun(@str2num,table2array(camels_CL_data(45,2:end))');
aridity_tmpa = cellfun(@str2num,table2array(camels_CL_data(46,2:end))');
p_seasonality_cr2met = cellfun(@str2num,table2array(camels_CL_data(47,2:end))');
p_seasonality_chirps = cellfun(@str2num,table2array(camels_CL_data(48,2:end))');
p_seasonality_mswep = cellfun(@str2num,table2array(camels_CL_data(49,2:end))');
p_seasonality_tmpa = cellfun(@str2num,table2array(camels_CL_data(50,2:end))');
frac_snow_cr2met = cellfun(@str2num,table2array(camels_CL_data(51,2:end))');
frac_snow_chirps = cellfun(@str2num,table2array(camels_CL_data(52,2:end))');
frac_snow_mswep = cellfun(@str2num,table2array(camels_CL_data(53,2:end))');
frac_snow_tmpa = cellfun(@str2num,table2array(camels_CL_data(54,2:end))');
high_prec_freq_cr2met = cellfun(@str2num,table2array(camels_CL_data(55,2:end))');
high_prec_freq_chirps = cellfun(@str2num,table2array(camels_CL_data(56,2:end))');
high_prec_freq_mswep = cellfun(@str2num,table2array(camels_CL_data(57,2:end))');
high_prec_freq_tmpa = cellfun(@str2num,table2array(camels_CL_data(58,2:end))');
high_prec_dur_cr2met = cellfun(@str2num,table2array(camels_CL_data(59,2:end))');
high_prec_dur_chirps = cellfun(@str2num,table2array(camels_CL_data(60,2:end))');
high_prec_dur_mswep = cellfun(@str2num,table2array(camels_CL_data(61,2:end))');
high_prec_dur_tmpa = cellfun(@str2num,table2array(camels_CL_data(62,2:end))');
high_prec_timing_cr2met = table2array(camels_CL_data(63,2:end))';
high_prec_timing_chirps = table2array(camels_CL_data(64,2:end))';
high_prec_timing_mswep = table2array(camels_CL_data(65,2:end))';
high_prec_timing_tmpa = table2array(camels_CL_data(66,2:end))';
low_prec_freq_cr2met = cellfun(@str2num,table2array(camels_CL_data(67,2:end))');
low_prec_freq_chirps = cellfun(@str2num,table2array(camels_CL_data(68,2:end))');
low_prec_freq_mswep = cellfun(@str2num,table2array(camels_CL_data(69,2:end))');
low_prec_freq_tmpa = cellfun(@str2num,table2array(camels_CL_data(70,2:end))');
low_prec_dur_cr2met = cellfun(@str2num,table2array(camels_CL_data(71,2:end))');
low_prec_dur_chirps = cellfun(@str2num,table2array(camels_CL_data(72,2:end))');
low_prec_dur_mswep = cellfun(@str2num,table2array(camels_CL_data(73,2:end))');
low_prec_dur_tmpa = cellfun(@str2num,table2array(camels_CL_data(74,2:end))');
low_prec_timing_cr2met = table2array(camels_CL_data(75,2:end))';
low_prec_timing_chirps = table2array(camels_CL_data(76,2:end))';
low_prec_timing_mswep = table2array(camels_CL_data(77,2:end))';
low_prec_timing_tmpa = table2array(camels_CL_data(78,2:end))';
p_mean_spread = cellfun(@str2num,table2array(camels_CL_data(79,2:end))');
% hydrological signatures
q_mean = cellfun(@str2num,table2array(camels_CL_data(80,2:end))');
runoff_ratio_cr2met = cellfun(@str2num,table2array(camels_CL_data(81,2:end))');
runoff_ratio_chirps = cellfun(@str2num,table2array(camels_CL_data(82,2:end))');
runoff_ratio_mswep = cellfun(@str2num,table2array(camels_CL_data(82,2:end))');
runoff_ratio_tmpa = cellfun(@str2num,table2array(camels_CL_data(83,2:end))');
stream_elas_cr2met = cellfun(@str2num,table2array(camels_CL_data(84,2:end))');
stream_elas_chirps = cellfun(@str2num,table2array(camels_CL_data(85,2:end))');
stream_elas_mswep = cellfun(@str2num,table2array(camels_CL_data(86,2:end))');
stream_elas_tmpa = cellfun(@str2num,table2array(camels_CL_data(87,2:end))');
slope_fdc = cellfun(@str2num,table2array(camels_CL_data(88,2:end))');
baseflow_index = cellfun(@str2num,table2array(camels_CL_data(89,2:end))');
hfd_mean = cellfun(@str2num,table2array(camels_CL_data(90,2:end))');
q95 = cellfun(@str2num,table2array(camels_CL_data(91,2:end))');
q5 = cellfun(@str2num,table2array(camels_CL_data(92,2:end))');
high_q_freq = cellfun(@str2num,table2array(camels_CL_data(93,2:end))');
high_q_dur = cellfun(@str2num,table2array(camels_CL_data(94,2:end))');
low_q_freq = cellfun(@str2num,table2array(camels_CL_data(95,2:end))');
low_q_dur = cellfun(@str2num,table2array(camels_CL_data(96,2:end))');
zero_q_freq = cellfun(@str2num,table2array(camels_CL_data(97,2:end))');
swe_ratio = cellfun(@str2num,table2array(camels_CL_data(98,2:end))');
% human intervention
sur_rights_n = cellfun(@str2num,table2array(camels_CL_data(99,2:end))');
sur_rights_flow = cellfun(@str2num,table2array(camels_CL_data(100,2:end))');
interv_degree = cellfun(@str2num,table2array(camels_CL_data(101,2:end))');
gw_rights_n = cellfun(@str2num,table2array(camels_CL_data(102,2:end))');
gw_rights_flow = cellfun(@str2num,table2array(camels_CL_data(103,2:end))');
big_dam = cellfun(@str2num,table2array(camels_CL_data(104,2:end))');

%% load time series
flow_perc_complete = NaN(length(gauge_id),1);
P = cell(length(gauge_id),1); % precipitation
PET = cell(length(gauge_id),1); % potential evapotranspiration
Q = cell(length(gauge_id),1); % streamflow
T = cell(length(gauge_id),1); % temperature

Q_table = readtable(...
    strcat(path_time_series,'\3_CAMELScl_streamflow_mm\3_CAMELScl_streamflow_mm.txt'));
Q_date = datenum(table2array(Q_table(2:end,1)));
% replace empty entries with NaN
for i = 1:size(Q_table,2)
    char1 = ' ';
    char2 = 'NaN';
    eval(['Q_table.Var',num2str(i),...
        '(strcmp(Q_table.Var',num2str(i),',char1)) = {char2};'])
end

P_table = readtable(...
    strcat(path_time_series,'\6_CAMELScl_precip_mswep\6_CAMELScl_precip_mswep.txt'));
P_date = datenum(table2array(P_table(2:end,1)));
% replace empty entries with NaN
for i = 1:size(P_table,2)
    char1 = ' ';
    char2 = 'NaN';
    eval(['P_table.Var',num2str(i),...
        '(strcmp(P_table.Var',num2str(i),',char1)) = {char2};'])
end

PET_table = readtable(...
    strcat(path_time_series,'\12_CAMELScl_pet_hargreaves\12_CAMELScl_pet_hargreaves.txt'));
PET_date = datenum(table2array(PET_table(2:end,1)));
% replace empty entries with NaN
for i = 1:size(PET_table,2)
    char1 = ' ';
    char2 = 'NaN';
    eval(['PET_table.Var',num2str(i),...
        '(strcmp(PET_table.Var',num2str(i),',char1)) = {char2};'])
end

T_table = readtable(...
    strcat(path_time_series,'\10_CAMELScl_tmean_cr2met\10_CAMELScl_tmean_cr2met.txt'));
T_date = datenum(table2array(T_table(2:end,1)));
% replace empty entries with NaN
for i = 1:size(T_table,2)
    char1 = ' ';
    char2 = 'NaN';
    eval(['T_table.Var',num2str(i),...
        '(strcmp(T_table.Var',num2str(i),',char1)) = {char2};'])
end

% loop over all catchments
for i = 1:length(gauge_id)
    
    if mod(i,10) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,length(gauge_id))
    end
    
    %ID = gauge_id(i);
    [P{i}, PET{i}, Q{i}, T{i}] = loadCatchmentCAMELSCL(...
        i,Q_table,Q_date,P_table,P_date,PET_table,PET_date,T_table,T_date);
    flow_perc_complete(i) = 100*(1-sum(isnan(Q{i}(:,2)))./length(Q{i}(:,2)));
    
end

%% create structure
% topography and metadata
CAMELS_CL_data.gauge_id = gauge_id;
CAMELS_CL_data.gauge_name = gauge_name;
CAMELS_CL_data.gauge_lat = gauge_lat;
CAMELS_CL_data.gauge_lon = gauge_lon;
% CAMELS_CL_data.record_period_start = record_period_start;
% CAMELS_CL_data.record_period_end = record_period_end;
% CAMELS_CL_data.n_obs = n_obs;
CAMELS_CL_data.area = area;
CAMELS_CL_data.elev_gauge = elev_gauge;
CAMELS_CL_data.elev_mean = elev_mean;
CAMELS_CL_data.elev_med = elev_med;
CAMELS_CL_data.elev_max = elev_max;
CAMELS_CL_data.elev_min = elev_min;
CAMELS_CL_data.slope_mean = slope_mean;
CAMELS_CL_data.nested_inner = nested_inner;
CAMELS_CL_data.nested_outer = nested_outer;
CAMELS_CL_data.location_type = location_type;
% geology
CAMELS_CL_data.geol_class_1st = geol_class_1st;
CAMELS_CL_data.geol_class_1st_frac = geol_class_1st_frac;
CAMELS_CL_data.geol_class_2nd = geol_class_2nd;
CAMELS_CL_data.geol_class_2nd_frac = geol_class_2nd_frac;
CAMELS_CL_data.carb_rocks_frac = carb_rocks_frac;
% land use
CAMELS_CL_data.crop_frac = crop_frac;
CAMELS_CL_data.nf_frac = nf_frac;
CAMELS_CL_data.fp_frac = fp_frac;
CAMELS_CL_data.grass_frac = grass_frac;
CAMELS_CL_data.shrub_frac = shrub_frac;
CAMELS_CL_data.wet_frac = wet_frac;
CAMELS_CL_data.imp_frac = imp_frac;
CAMELS_CL_data.lc_barren = lc_barren;
CAMELS_CL_data.snow_frac = snow_frac;
CAMELS_CL_data.lc_glacier = lc_glacier;
CAMELS_CL_data.fp_nf_index = fp_nf_index;
CAMELS_CL_data.forest_frac = forest_frac;
CAMELS_CL_data.dom_land_cover = dom_land_cover;
CAMELS_CL_data.dom_land_cover_frac = dom_land_cover_frac;
CAMELS_CL_data.land_cover_missing = land_cover_missing;
% climate
CAMELS_CL_data.p_mean_cr2met = p_mean_cr2met;
CAMELS_CL_data.p_mean_chirps = p_mean_chirps;
CAMELS_CL_data.p_mean_mswep = p_mean_mswep;
CAMELS_CL_data.p_mean_tmpa = p_mean_tmpa;
CAMELS_CL_data.pet_mean = pet_mean;
CAMELS_CL_data.aridity_cr2met = aridity_cr2met;
CAMELS_CL_data.aridity_chirps = aridity_chirps;
CAMELS_CL_data.aridity_mswep = aridity_mswep;
CAMELS_CL_data.aridity_tmpa = aridity_tmpa;
CAMELS_CL_data.p_seasonality_cr2met = p_seasonality_cr2met;
CAMELS_CL_data.p_seasonality_chirps = p_seasonality_chirps;
CAMELS_CL_data.p_seasonality_mswep = p_seasonality_mswep;
CAMELS_CL_data.p_seasonality_tmpa = p_seasonality_tmpa;
CAMELS_CL_data.frac_snow_cr2met = frac_snow_cr2met;
CAMELS_CL_data.frac_snow_chirps = frac_snow_chirps;
CAMELS_CL_data.frac_snow_mswep = frac_snow_mswep;
CAMELS_CL_data.frac_snow_tmpa = frac_snow_tmpa;
CAMELS_CL_data.high_prec_freq_cr2met = high_prec_freq_cr2met;
CAMELS_CL_data.high_prec_freq_chirps = high_prec_freq_chirps;
CAMELS_CL_data.high_prec_freq_mswep = high_prec_freq_mswep;
CAMELS_CL_data.high_prec_freq_tmpa = high_prec_freq_tmpa;
CAMELS_CL_data.high_prec_dur_cr2met = high_prec_dur_cr2met;
CAMELS_CL_data.high_prec_dur_chirps = high_prec_dur_chirps;
CAMELS_CL_data.high_prec_dur_mswep = high_prec_dur_mswep;
CAMELS_CL_data.high_prec_dur_tmpa = high_prec_dur_tmpa;
CAMELS_CL_data.high_prec_timing_cr2met = high_prec_timing_cr2met;
CAMELS_CL_data.high_prec_timing_chirps = high_prec_timing_chirps;
CAMELS_CL_data.high_prec_timing_mswep = high_prec_timing_mswep;
CAMELS_CL_data.high_prec_timing_tmpa = high_prec_timing_tmpa;
CAMELS_CL_data.low_prec_freq_cr2met = low_prec_freq_cr2met;
CAMELS_CL_data.low_prec_freq_chirps = low_prec_freq_chirps;
CAMELS_CL_data.low_prec_freq_mswep = low_prec_freq_mswep;
CAMELS_CL_data.low_prec_freq_tmpa = low_prec_freq_tmpa;
CAMELS_CL_data.low_prec_dur_cr2met = low_prec_dur_cr2met;
CAMELS_CL_data.low_prec_dur_chirps = low_prec_dur_chirps;
CAMELS_CL_data.low_prec_dur_mswep = low_prec_dur_mswep;
CAMELS_CL_data.low_prec_dur_tmpa = low_prec_dur_tmpa;
CAMELS_CL_data.low_prec_timing_cr2met = low_prec_timing_cr2met;
CAMELS_CL_data.low_prec_timing_chirps = low_prec_timing_chirps;
CAMELS_CL_data.low_prec_timing_mswep = low_prec_timing_mswep;
CAMELS_CL_data.low_prec_timing_tmpa = low_prec_timing_tmpa;
CAMELS_CL_data.p_mean_spread = p_mean_spread;
% hydrological signatures
CAMELS_CL_data.q_mean = q_mean;
CAMELS_CL_data.runoff_ratio_cr2met = runoff_ratio_cr2met;
CAMELS_CL_data.runoff_ratio_chirps = runoff_ratio_chirps;
CAMELS_CL_data.runoff_ratio_mswep =runoff_ratio_mswep ;
CAMELS_CL_data.runoff_ratio_tmpa = runoff_ratio_tmpa;
CAMELS_CL_data.stream_elas_cr2met = stream_elas_cr2met;
CAMELS_CL_data.stream_elas_chirps = stream_elas_chirps;
CAMELS_CL_data.stream_elas_mswep = stream_elas_mswep;
CAMELS_CL_data.stream_elas_tmpa = stream_elas_tmpa;
CAMELS_CL_data.slope_fdc = slope_fdc;
CAMELS_CL_data.baseflow_index = baseflow_index;
CAMELS_CL_data.hfd_mean = hfd_mean;
CAMELS_CL_data.q95 = q95;
CAMELS_CL_data.q5 = q5;
CAMELS_CL_data.high_q_freq = high_q_freq;
CAMELS_CL_data.high_q_dur = high_q_dur;
CAMELS_CL_data.low_q_freq = low_q_freq;
CAMELS_CL_data.low_q_dur = low_q_dur;
CAMELS_CL_data.zero_q_freq = zero_q_freq;
CAMELS_CL_data.swe_ratio = swe_ratio;
% human intervention
CAMELS_CL_data.sur_rights_n = sur_rights_n;
CAMELS_CL_data.sur_rights_flow = sur_rights_flow;
CAMELS_CL_data.interv_degree = interv_degree;
CAMELS_CL_data.gw_rights_n = gw_rights_n;
CAMELS_CL_data.gw_rights_flow = gw_rights_flow;
CAMELS_CL_data.big_dam = big_dam;

% hydro-meteorological time series
CAMELS_CL_data.flow_perc_complete = flow_perc_complete;
CAMELS_CL_data.P = P;
CAMELS_CL_data.PET = PET;
CAMELS_CL_data.Q = Q;
CAMELS_CL_data.T = T;

% save file to table
save('./CAMELS_Matlab/Data/CAMELS_CL_data.mat','CAMELS_CL_data')

