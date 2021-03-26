function [CAMELS_CL_data] = saveCAMELSstruct_CL(save_struct)
%saveCAMELSstruct_CL Creates struct file with CAMELS CL data.
%   - Loads hydro-meteorological time series and catchment attributes
%   - Timeseries are loaded for the period in which all data are available 
%   - Uses local paths (which contain large CAMELS files)
%   - Uses MSWEP P and HARGREAVES PET data for time series
%   - Data can be found at: https://doi.pangaea.de/10.1594/PANGAEA.894885?format=html#download
%
%   INPUT
%   save_struct: whether to save struct file or not
%
%   OUTPUT
%   CAMELS_CL_data: struct file with CAMELS CL data
%
%   References
%   Alvarez-Garreton, C., Mendoza, P.A., Boisier, J.P., Addor, N., 
%   Galleguillos, M., Zambrano-Bigiarini, M., Lara, A., Puelma, C., Cortes, 
%   G., Garreaud, R. and McPhee, J., 2018. The CAMELS-CL dataset: catchment 
%   attributes and meteorology for large sample studies–Chile dataset. 
%   Hydrology and Earth System Sciences, 22(11), pp.5817-5846.
%
%   Copyright (C) 2021
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 1
    save_struct = false;
end

%% Specify paths
% We have to be above the CAMELS_Matlab directory and CAMELS data should be 
% stored in a folder named CAMELS_CL. The following folders are required 
% (each contain one txt file):
% CAMELS_CL/1_CAMELScl_attributes
% CAMELS_CL/3_CAMELScl_streamflow_mm
% CAMELS_CL/5_CAMELScl_precip_chirps
% CAMELS_CL/6_CAMELScl_precip_mswep
% CAMELS_CL/10_CAMELScl_tmean_cr2met
% CAMELS_CL/11_CAMELScl_pet_8d_modis
% CAMELS_CL/12_CAMELScl_pet_hargreaves

path_catchment_attributes = "CAMELS_CL/1_CAMELScl_attributes/";
path_time_series = "CAMELS_CL/";

if ~(exist(path_catchment_attributes) == 7)
    error('Cannot find local path. You can download CAMELS CL from https://doi.pangaea.de/10.1594/PANGAEA.894885?format=html#download.')
elseif ~(exist(path_time_series) == 7)
    error('Cannot find local path. You can download CAMELS CL from https://doi.pangaea.de/10.1594/PANGAEA.894885?format=html#download.')
end

%% Load catchment attributes
% We first load the catchment attribute data which are saved in several csv
% files.

% Note: I could not find a better way to load these data (vertical headers
% are difficult to handle...). Please let me know if there's a more elegant
% way to load these data.
fmt = repmat('%s',1,517);
fileID = fopen(strcat(path_catchment_attributes,'1_CAMELScl_attributes.txt'),'r');
file = fread(fileID,'*char');
file = strrep(file','NA','NaN'); % replace NA with NaN
file = strrep(file,'"',' '); % replace " with space
file = strrep(file,'\',' '); % replace \ with space
fclose(fileID);
fileID = fopen(strcat(path_catchment_attributes,'data_tmp.txt'),'w');
fprintf(fileID,file);
fclose(fileID);

% all attributes
camels_CL_attribute_data_str = readtable(strcat(path_catchment_attributes,'data_tmp.txt'),'Format',fmt);
camels_CL_attribute_data = readtable(strcat(path_catchment_attributes,'1_CAMELScl_attributes.txt'));

% todo: delete temporary file
% delete(strcat(path_catchment_attributes,'data_tmp.txt'));

% replace NA with NaN 
% for i = 1:size(camels_CL_data,2)
%     char1 = 'NA';
%     char2 = 'NaN';
%     eval(['camels_CL_data.Var',num2str(i),...
%         '(strcmp(camels_CL_data.Var',num2str(i),',char1)) = {char2};'])
% end

% We now add the catchment attributes and metadata to the struct file.

% topography and metadata
CAMELS_CL_data.gauge_id = table2array(camels_CL_attribute_data(1,2:end))';
CAMELS_CL_data.gauge_name = table2array(camels_CL_attribute_data_str(2,2:end))';
CAMELS_CL_data.gauge_lat = table2array(camels_CL_attribute_data(3,2:end))';
CAMELS_CL_data.gauge_lon = table2array(camels_CL_attribute_data(4,2:end))';
% record_period_start
% record_period_end
% n_obs
CAMELS_CL_data.area = table2array(camels_CL_attribute_data(8,2:end))';
CAMELS_CL_data.elev_gauge = table2array(camels_CL_attribute_data(9,2:end))';
CAMELS_CL_data.elev_mean = table2array(camels_CL_attribute_data(10,2:end))';
CAMELS_CL_data.elev_med = table2array(camels_CL_attribute_data(11,2:end))';
CAMELS_CL_data.elev_max = table2array(camels_CL_attribute_data(12,2:end))';
CAMELS_CL_data.elev_min = table2array(camels_CL_attribute_data(13,2:end))';
CAMELS_CL_data.slope_mean = table2array(camels_CL_attribute_data(14,2:end))';
CAMELS_CL_data.nested_inner = table2array(camels_CL_attribute_data(15,2:end))';
CAMELS_CL_data.nested_outer = table2array(camels_CL_attribute_data(16,2:end))';
CAMELS_CL_data.location_type = table2array(camels_CL_attribute_data_str(17,2:end))';

% geology
CAMELS_CL_data.geol_class_1st = table2array(camels_CL_attribute_data(18,2:end))';
CAMELS_CL_data.geol_class_1st_frac = table2array(camels_CL_attribute_data(19,2:end))';
CAMELS_CL_data.geol_class_2nd = table2array(camels_CL_attribute_data(1,20:end))';
CAMELS_CL_data.geol_class_2nd_frac = table2array(camels_CL_attribute_data(21,2:end))';
CAMELS_CL_data.carb_rocks_frac = table2array(camels_CL_attribute_data(22,2:end))';

% land use
CAMELS_CL_data.crop_frac = table2array(camels_CL_attribute_data(23,2:end))';
CAMELS_CL_data.nf_frac = table2array(camels_CL_attribute_data(24,2:end))';
CAMELS_CL_data.fp_frac = table2array(camels_CL_attribute_data(25,2:end))';
CAMELS_CL_data.grass_frac = table2array(camels_CL_attribute_data(26,2:end))';
CAMELS_CL_data.shrub_frac = table2array(camels_CL_attribute_data(27,2:end))';
CAMELS_CL_data.wet_frac = table2array(camels_CL_attribute_data(28,2:end))';
CAMELS_CL_data.imp_frac = table2array(camels_CL_attribute_data(29,2:end))';
CAMELS_CL_data.lc_barren = table2array(camels_CL_attribute_data(30,2:end))';
CAMELS_CL_data.snow_frac = table2array(camels_CL_attribute_data(31,2:end))';
CAMELS_CL_data.lc_glacier = table2array(camels_CL_attribute_data(32,2:end))';
CAMELS_CL_data.fp_nf_index = table2array(camels_CL_attribute_data(33,2:end))';
CAMELS_CL_data.forest_frac = table2array(camels_CL_attribute_data(34,2:end))';
CAMELS_CL_data.dom_land_cover = table2array(camels_CL_attribute_data(35,2:end))';
CAMELS_CL_data.dom_land_cover_frac = table2array(camels_CL_attribute_data(36,2:end))';
CAMELS_CL_data.land_cover_missing = table2array(camels_CL_attribute_data(37,2:end))';

% climate
CAMELS_CL_data.p_mean_cr2met = table2array(camels_CL_attribute_data(38,2:end))';
CAMELS_CL_data.p_mean_chirps = table2array(camels_CL_attribute_data(39,2:end))';
CAMELS_CL_data.p_mean_mswep = table2array(camels_CL_attribute_data(40,2:end))';
CAMELS_CL_data.p_mean_tmpa = table2array(camels_CL_attribute_data(41,2:end))';
CAMELS_CL_data.pet_mean = table2array(camels_CL_attribute_data(42,2:end))';
CAMELS_CL_data.aridity_cr2met = table2array(camels_CL_attribute_data(43,2:end))';
CAMELS_CL_data.aridity_chirps = table2array(camels_CL_attribute_data(44,2:end))';
CAMELS_CL_data.aridity_mswep = table2array(camels_CL_attribute_data(45,2:end))';
CAMELS_CL_data.aridity_tmpa = table2array(camels_CL_attribute_data(46,2:end))';
CAMELS_CL_data.p_seasonality_cr2met = table2array(camels_CL_attribute_data(47,2:end))';
CAMELS_CL_data.p_seasonality_chirps = table2array(camels_CL_attribute_data(48,2:end))';
CAMELS_CL_data.p_seasonality_mswep = table2array(camels_CL_attribute_data(49,2:end))';
CAMELS_CL_data.p_seasonality_tmpa = table2array(camels_CL_attribute_data(50,2:end))';
CAMELS_CL_data.frac_snow_cr2met = table2array(camels_CL_attribute_data(51,2:end))';
CAMELS_CL_data.frac_snow_chirps = table2array(camels_CL_attribute_data(52,2:end))';
CAMELS_CL_data.frac_snow_mswep = table2array(camels_CL_attribute_data(53,2:end))';
CAMELS_CL_data.frac_snow_tmpa = table2array(camels_CL_attribute_data(54,2:end))';
CAMELS_CL_data.high_prec_freq_cr2met = table2array(camels_CL_attribute_data(55,2:end))';
CAMELS_CL_data.high_prec_freq_chirps = table2array(camels_CL_attribute_data(56,2:end))';
CAMELS_CL_data.high_prec_freq_mswep = table2array(camels_CL_attribute_data(57,2:end))';
CAMELS_CL_data.high_prec_freq_tmpa = table2array(camels_CL_attribute_data(58,2:end))';
CAMELS_CL_data.high_prec_dur_cr2met = table2array(camels_CL_attribute_data(59,2:end))';
CAMELS_CL_data.high_prec_dur_chirps = table2array(camels_CL_attribute_data(60,2:end))';
CAMELS_CL_data.high_prec_dur_mswep = table2array(camels_CL_attribute_data(61,2:end))';
CAMELS_CL_data.high_prec_dur_tmpa = table2array(camels_CL_attribute_data(62,2:end))';
CAMELS_CL_data.high_prec_timing_cr2met = table2array(camels_CL_attribute_data(63,2:end))';
CAMELS_CL_data.high_prec_timing_chirps = table2array(camels_CL_attribute_data(64,2:end))';
CAMELS_CL_data.high_prec_timing_mswep = table2array(camels_CL_attribute_data(65,2:end))';
CAMELS_CL_data.high_prec_timing_tmpa = table2array(camels_CL_attribute_data(66,2:end))';
CAMELS_CL_data.low_prec_freq_cr2met = table2array(camels_CL_attribute_data(67,2:end))';
CAMELS_CL_data.low_prec_freq_chirps = table2array(camels_CL_attribute_data(68,2:end))';
CAMELS_CL_data.low_prec_freq_mswep = table2array(camels_CL_attribute_data(69,2:end))';
CAMELS_CL_data.low_prec_freq_tmpa = table2array(camels_CL_attribute_data(70,2:end))';
CAMELS_CL_data.low_prec_dur_cr2met = table2array(camels_CL_attribute_data(71,2:end))';
CAMELS_CL_data.low_prec_dur_chirps = table2array(camels_CL_attribute_data(72,2:end))';
CAMELS_CL_data.low_prec_dur_mswep = table2array(camels_CL_attribute_data(73,2:end))';
CAMELS_CL_data.low_prec_dur_tmpa = table2array(camels_CL_attribute_data(74,2:end))';
CAMELS_CL_data.low_prec_timing_cr2met = table2array(camels_CL_attribute_data(75,2:end))';
CAMELS_CL_data.low_prec_timing_chirps = table2array(camels_CL_attribute_data(76,2:end))';
CAMELS_CL_data.low_prec_timing_mswep = table2array(camels_CL_attribute_data(77,2:end))';
CAMELS_CL_data.low_prec_timing_tmpa = table2array(camels_CL_attribute_data(78,2:end))';
CAMELS_CL_data.p_mean_spread = table2array(camels_CL_attribute_data(79,2:end))';

% hydrological signatures
CAMELS_CL_data.q_mean = table2array(camels_CL_attribute_data(80,2:end))';
CAMELS_CL_data.runoff_ratio_cr2met = table2array(camels_CL_attribute_data(81,2:end))';
CAMELS_CL_data.runoff_ratio_chirps = table2array(camels_CL_attribute_data(82,2:end))';
CAMELS_CL_data.runoff_ratio_mswep = table2array(camels_CL_attribute_data(83,2:end))';
CAMELS_CL_data.runoff_ratio_tmpa = table2array(camels_CL_attribute_data(84,2:end))';
CAMELS_CL_data.stream_elas_cr2met = table2array(camels_CL_attribute_data(85,2:end))';
CAMELS_CL_data.stream_elas_chirps = table2array(camels_CL_attribute_data(86,2:end))';
CAMELS_CL_data.stream_elas_mswep = table2array(camels_CL_attribute_data(87,2:end))';
CAMELS_CL_data.stream_elas_tmpa = table2array(camels_CL_attribute_data(88,2:end))';
CAMELS_CL_data.slope_fdc = table2array(camels_CL_attribute_data(89,2:end))';
CAMELS_CL_data.baseflow_index = table2array(camels_CL_attribute_data(90,2:end))';
CAMELS_CL_data.hfd_mean = table2array(camels_CL_attribute_data(91,2:end))';
CAMELS_CL_data.q95 = table2array(camels_CL_attribute_data(92,2:end))';
CAMELS_CL_data.q5 = table2array(camels_CL_attribute_data(93,2:end))';
CAMELS_CL_data.high_q_freq = table2array(camels_CL_attribute_data(94,2:end))';
CAMELS_CL_data.high_q_dur = table2array(camels_CL_attribute_data(95,2:end))';
CAMELS_CL_data.low_q_freq = table2array(camels_CL_attribute_data(96,2:end))';
CAMELS_CL_data.low_q_dur = table2array(camels_CL_attribute_data(97,2:end))';
CAMELS_CL_data.zero_q_freq = table2array(camels_CL_attribute_data(98,2:end))';
CAMELS_CL_data.swe_ratio = table2array(camels_CL_attribute_data(99,2:end))';

% human intervention
CAMELS_CL_data.sur_rights_n = table2array(camels_CL_attribute_data(100,2:end))';
CAMELS_CL_data.sur_rights_flow = table2array(camels_CL_attribute_data(101,2:end))';
CAMELS_CL_data.interv_degree = table2array(camels_CL_attribute_data(102,2:end))';
CAMELS_CL_data.gw_rights_n = table2array(camels_CL_attribute_data(103,2:end))';
CAMELS_CL_data.gw_rights_flow = table2array(camels_CL_attribute_data(104,2:end))';
CAMELS_CL_data.big_dam = table2array(camels_CL_attribute_data(105,2:end))';

%% Load hydro-meteorological time series
% To extract the time series, we loop over all catchments. We also
% calculate the completeness of the flow records.
flow_perc_complete = NaN(length(CAMELS_CL_data.gauge_id),1);
P = cell(length(CAMELS_CL_data.gauge_id),1); % precipitation
PET = cell(length(CAMELS_CL_data.gauge_id),1); % potential evapotranspiration
Q = cell(length(CAMELS_CL_data.gauge_id),1); % streamflow
T = cell(length(CAMELS_CL_data.gauge_id),1); % temperature

Q_table = readtable(...
    strcat(path_time_series,'/3_CAMELScl_streamflow_mm/3_CAMELScl_streamflow_mm.txt'));
Q_date = datenum(table2array(Q_table(2:end,1)));
% replace empty entries with NaN
% for i = 1:size(Q_table,2)
%     char1 = ' ';
%     char2 = 'NaN';
%     eval(['Q_table.Var',num2str(i),...
%         '(strcmp(Q_table.Var',num2str(i),',char1)) = {char2};'])
% end

P_table = readtable(...
    strcat(path_time_series,'/6_CAMELScl_precip_mswep/6_CAMELScl_precip_mswep.txt'));
P_date = datenum(table2array(P_table(2:end,1)));
% replace empty entries with NaN
% for i = 1:size(P_table,2)
%     char1 = ' ';
%     char2 = 'NaN';
%     eval(['P_table.Var',num2str(i),...
%         '(strcmp(P_table.Var',num2str(i),',char1)) = {char2};'])
% end

PET_table = readtable(...
    strcat(path_time_series,'/12_CAMELScl_pet_hargreaves/12_CAMELScl_pet_hargreaves.txt'));
PET_date = datenum(table2array(PET_table(2:end,1)));
% replace empty entries with NaN
% for i = 1:size(PET_table,2)
%     char1 = ' ';
%     char2 = 'NaN';
%     eval(['PET_table.Var',num2str(i),...
%         '(strcmp(PET_table.Var',num2str(i),',char1)) = {char2};'])
% end

T_table = readtable(...
    strcat(path_time_series,'/10_CAMELScl_tmean_cr2met/10_CAMELScl_tmean_cr2met.txt'));
T_date = datenum(table2array(T_table(2:end,1)));
% replace empty entries with NaN
% for i = 1:size(T_table,2)
%     char1 = ' ';
%     char2 = 'NaN';
%     eval(['T_table.Var',num2str(i),...
%         '(strcmp(T_table.Var',num2str(i),',char1)) = {char2};'])
% end

fprintf('Loading catchment data (CL)...\n') 
for i = 1:length(CAMELS_CL_data.gauge_id) % loop over all catchments
    
    if mod(i,100) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,length(CAMELS_CL_data.gauge_id))
    end
    
    [P{i}, PET{i}, Q{i}, T{i}] = loadCatchmentCAMELS_CL(...
        i,Q_table,Q_date,P_table,P_date,PET_table,PET_date,T_table,T_date);
    flow_perc_complete(i) = 100*(1-sum(isnan(Q{i}(:,2)))./length(Q{i}(:,2)));
    
end

% add hydro-meteorological time series to struct file
CAMELS_CL_data.flow_perc_complete = flow_perc_complete;
CAMELS_CL_data.P = P;
CAMELS_CL_data.PET = PET;
CAMELS_CL_data.Q = Q;
CAMELS_CL_data.T = T;

% save the struct file
if save_struct
    save('CAMELS_Matlab/Data/CAMELS_CL_data.mat','-struct','CAMELS_CL_data')
end

end
