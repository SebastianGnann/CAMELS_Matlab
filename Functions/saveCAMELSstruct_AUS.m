function [CAMELS_AUS_data] = saveCAMELSstruct_AUS(save_struct)
%saveCAMELSstruct_AUS Creates struct file with CAMELS AUS data.
%   - Loads hydro-meteorological time series and catchment attributes
%   - Timeseries are loaded for the period in which all data are available 
%   - Uses local paths (which contain large CAMELS files)
%   - Loads AWAP precipitation, Morton SILO PET, and AWAP temperature
%   - Data can be found at: https://doi.pangaea.de/10.1594/PANGAEA.921850?format=html#download
%   - Data paper is still under review
%
%   INPUT
%   save_struct: whether to save struct file or not
%
%   OUTPUT
%   CAMELS_AUS_data: struct file with CAMELS AUS data
%
%   References
%   Fowler, K.J., Acharya, S.C., Addor, N., Chou, C. and Peel, M.C., 2021. 
%   CAMELS-AUS: Hydrometeorological time series and landscape attributes 
%   for 222 catchments in Australia. Earth System Science Data Discussions, 
%   pp.1-30.
%
%   Copyright (C) 2021
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 1
    save_struct = false;
end


%% Specify paths
% We have to be above the CAMELS_Matlab directory and CAMELS data should be 
% stored in a folder named CAMELS_AUS. The following files/folders are  
% required:
% CAMELS_AUS/03_streamflow (contains streamflow time series)
% CAMELS_AUS/05_hydrometeorology (contains meteorological time series)
% CAMELS_AUS_Attributes-Indices_MasterTable (contains catchment attributes)

path_catchment_attributes = strcat("CAMELS_AUS/CAMELS_AUS_Attributes-Indices_MasterTable.csv");
path = "CAMELS_AUS/"; 

if ~(exist(path_catchment_attributes) == 2)
    error('Cannot find local path. You can download CAMELS AUS from https://doi.pangaea.de/10.1594/PANGAEA.921850?format=html#download.')
elseif ~(exist(path) == 7)
    error('Cannot find local path. You can download CAMELS AUS from https://doi.pangaea.de/10.1594/PANGAEA.921850?format=html#download.')
end

%% Load catchment attributes
% We first load the catchment attribute data which are saved in several txt 
% files.

% all attributes
opts = detectImportOptions(path_catchment_attributes);
opts = setvartype(opts, "station_id", 'char');  
camels_AUS_attribute_data = readtable(path_catchment_attributes,opts);

% We now add the catchment attributes and metadata to the struct file.

% geography and metadata
CAMELS_AUS_data.station_id = camels_AUS_attribute_data.station_id;
gauge_id = string(camels_AUS_attribute_data.station_id); 
match = ["A","B","D","G"]; gauge_id = double(erase(gauge_id,match)); 
CAMELS_AUS_data.gauge_id = gauge_id; % make station IDs numeric
CAMELS_AUS_data.station_name = camels_AUS_attribute_data.station_name;
CAMELS_AUS_data.drainage_division = camels_AUS_attribute_data.drainage_division;
CAMELS_AUS_data.river_region = camels_AUS_attribute_data.river_region;
CAMELS_AUS_data.notes = camels_AUS_attribute_data.notes;
CAMELS_AUS_data.lat_outlet = camels_AUS_attribute_data.lat_outlet;
CAMELS_AUS_data.long_outlet = camels_AUS_attribute_data.long_outlet;
CAMELS_AUS_data.lat_centroid = camels_AUS_attribute_data.lat_centroid;
CAMELS_AUS_data.long_centroid = camels_AUS_attribute_data.long_centroid;
CAMELS_AUS_data.map_zone = camels_AUS_attribute_data.map_zone;
CAMELS_AUS_data.catchment_area = camels_AUS_attribute_data.catchment_area;
CAMELS_AUS_data.nested_status = camels_AUS_attribute_data.nested_status;
CAMELS_AUS_data.next_station_ds = camels_AUS_attribute_data.next_station_ds;
CAMELS_AUS_data.num_nested_within = camels_AUS_attribute_data.num_nested_within;
CAMELS_AUS_data.start_date = camels_AUS_attribute_data.start_date;
CAMELS_AUS_data.end_date = camels_AUS_attribute_data.end_date;
CAMELS_AUS_data.long_centroid = camels_AUS_attribute_data.long_centroid;
CAMELS_AUS_data.prop_missing_data = camels_AUS_attribute_data.prop_missing_data;

% streamflow uncertainty 								
CAMELS_AUS_data.q_uncert_num_curves = camels_AUS_attribute_data.q_uncert_num_curves;
CAMELS_AUS_data.q_uncert_n = camels_AUS_attribute_data.q_uncert_n;
CAMELS_AUS_data.q_uncert_q10 = camels_AUS_attribute_data.q_uncert_q10;
CAMELS_AUS_data.q_uncert_q10_upper = camels_AUS_attribute_data.q_uncert_q10_upper;
CAMELS_AUS_data.q_uncert_q10_lower = camels_AUS_attribute_data.q_uncert_q10_lower;
CAMELS_AUS_data.q_uncert_q50 = camels_AUS_attribute_data.q_uncert_q50;
CAMELS_AUS_data.q_uncert_q50_upper = camels_AUS_attribute_data.q_uncert_q50_upper;
CAMELS_AUS_data.q_uncert_q50_lower = camels_AUS_attribute_data.q_uncert_q50_lower;
CAMELS_AUS_data.q_uncert_q90 = camels_AUS_attribute_data.q_uncert_q90;
CAMELS_AUS_data.q_uncert_q90_upper = camels_AUS_attribute_data.q_uncert_q90_upper;
CAMELS_AUS_data.q_uncert_q90_lower = camels_AUS_attribute_data.q_uncert_q90_lower;

% climate indices						
CAMELS_AUS_data.p_mean = camels_AUS_attribute_data.p_mean;
CAMELS_AUS_data.pet_mean = camels_AUS_attribute_data.pet_mean;
CAMELS_AUS_data.aridity = camels_AUS_attribute_data.aridity;
CAMELS_AUS_data.p_seasonality = camels_AUS_attribute_data.p_seasonality;
CAMELS_AUS_data.frac_snow = camels_AUS_attribute_data.frac_snow;
CAMELS_AUS_data.high_prec_freq = camels_AUS_attribute_data.high_prec_freq;
CAMELS_AUS_data.high_prec_dur = camels_AUS_attribute_data.high_prec_dur;
CAMELS_AUS_data.high_prec_timing = camels_AUS_attribute_data.high_prec_timing;
CAMELS_AUS_data.low_prec_freq = camels_AUS_attribute_data.low_prec_freq;
CAMELS_AUS_data.low_prec_dur = camels_AUS_attribute_data.low_prec_dur;
CAMELS_AUS_data.low_prec_timing = camels_AUS_attribute_data.low_prec_timing;

% hydrological signatures											
CAMELS_AUS_data.q_mean = camels_AUS_attribute_data.q_mean;
CAMELS_AUS_data.runoff_ratio = camels_AUS_attribute_data.runoff_ratio;
CAMELS_AUS_data.stream_elas = camels_AUS_attribute_data.stream_elas;
CAMELS_AUS_data.slope_fdc = camels_AUS_attribute_data.slope_fdc;
CAMELS_AUS_data.baseflow_index = camels_AUS_attribute_data.baseflow_index;
CAMELS_AUS_data.hdf_mean = camels_AUS_attribute_data.hdf_mean;
CAMELS_AUS_data.Q5 = camels_AUS_attribute_data.Q5;
CAMELS_AUS_data.Q95 = camels_AUS_attribute_data.Q95;
CAMELS_AUS_data.high_q_freq = camels_AUS_attribute_data.high_q_freq;
CAMELS_AUS_data.high_q_dur = camels_AUS_attribute_data.high_q_dur;
CAMELS_AUS_data.low_q_freq = camels_AUS_attribute_data.low_q_freq;
CAMELS_AUS_data.low_q_dur = camels_AUS_attribute_data.zero_q_freq;

% geology and soils
CAMELS_AUS_data.geol_prim = camels_AUS_attribute_data.geol_prim;				
CAMELS_AUS_data.geol_prim_prop = camels_AUS_attribute_data.geol_prim_prop;				
CAMELS_AUS_data.geol_sec = camels_AUS_attribute_data.geol_sec;				
CAMELS_AUS_data.geol_sec_prop = camels_AUS_attribute_data.geol_sec_prop;				
CAMELS_AUS_data.unconsoldted = camels_AUS_attribute_data.unconsoldted;				
CAMELS_AUS_data.igneous = camels_AUS_attribute_data.igneous;			
CAMELS_AUS_data.silicsed = camels_AUS_attribute_data.silicsed;			
CAMELS_AUS_data.carbnatesed = camels_AUS_attribute_data.carbnatesed;			
CAMELS_AUS_data.othersed = camels_AUS_attribute_data.othersed;			
CAMELS_AUS_data.metamorph = camels_AUS_attribute_data.metamorph;			
CAMELS_AUS_data.sedvolc = camels_AUS_attribute_data.sedvolc;			
CAMELS_AUS_data.oldrock = camels_AUS_attribute_data.oldrock;		
CAMELS_AUS_data.claya = camels_AUS_attribute_data.claya;		
CAMELS_AUS_data.clayb = camels_AUS_attribute_data.clayb;		
CAMELS_AUS_data.sanda = camels_AUS_attribute_data.sanda;		
CAMELS_AUS_data.solum_thickness = camels_AUS_attribute_data.solum_thickness;
CAMELS_AUS_data.ksat = camels_AUS_attribute_data.ksat;
CAMELS_AUS_data.solpawhc = camels_AUS_attribute_data.solpawhc;

% topography and geometry
CAMELS_AUS_data.elev_min = camels_AUS_attribute_data.elev_min;				
CAMELS_AUS_data.elev_max = camels_AUS_attribute_data.elev_max;				
CAMELS_AUS_data.elev_mean = camels_AUS_attribute_data.elev_mean;				
CAMELS_AUS_data.elev_range = camels_AUS_attribute_data.elev_range;				
CAMELS_AUS_data.mean_slope_pct = camels_AUS_attribute_data.mean_slope_pct;				
CAMELS_AUS_data.upsdist = camels_AUS_attribute_data.upsdist;			
CAMELS_AUS_data.strdensity = camels_AUS_attribute_data.strdensity;		
CAMELS_AUS_data.strahler = camels_AUS_attribute_data.strahler;		
CAMELS_AUS_data.elongratio = camels_AUS_attribute_data.elongratio;		
CAMELS_AUS_data.relief = camels_AUS_attribute_data.relief;		
CAMELS_AUS_data.reliefratio = camels_AUS_attribute_data.reliefratio;		
CAMELS_AUS_data.mrvbf_prop_0 = camels_AUS_attribute_data.mrvbf_prop_0;		
CAMELS_AUS_data.mrvbf_prop_1 = camels_AUS_attribute_data.mrvbf_prop_1;		
CAMELS_AUS_data.mrvbf_prop_2 = camels_AUS_attribute_data.mrvbf_prop_2;	
CAMELS_AUS_data.mrvbf_prop_3 = camels_AUS_attribute_data.mrvbf_prop_3;	
CAMELS_AUS_data.mrvbf_prop_4 = camels_AUS_attribute_data.mrvbf_prop_4;	
CAMELS_AUS_data.mrvbf_prop_5 = camels_AUS_attribute_data.mrvbf_prop_5;	
CAMELS_AUS_data.mrvbf_prop_6 = camels_AUS_attribute_data.mrvbf_prop_6;	
CAMELS_AUS_data.mrvbf_prop_7 = camels_AUS_attribute_data.mrvbf_prop_7;	
CAMELS_AUS_data.mrvbf_prop_8 = camels_AUS_attribute_data.mrvbf_prop_8;	
CAMELS_AUS_data.mrvbf_prop_9 = camels_AUS_attribute_data.mrvbf_prop_9;
CAMELS_AUS_data.confinement = camels_AUS_attribute_data.confinement;

% land cover and vegetation
CAMELS_AUS_data.lc01_extracti = camels_AUS_attribute_data.lc01_extracti;
CAMELS_AUS_data.lc03_waterbo = camels_AUS_attribute_data.lc03_waterbo;
CAMELS_AUS_data.lc04_saltlak = camels_AUS_attribute_data.lc04_saltlak;
CAMELS_AUS_data.lc05_irrcrop = camels_AUS_attribute_data.lc05_irrcrop;
CAMELS_AUS_data.lc06_irrpast = camels_AUS_attribute_data.lc06_irrpast;
CAMELS_AUS_data.lc07_irrsuga = camels_AUS_attribute_data.lc07_irrsuga;
CAMELS_AUS_data.lc08_rfcropp = camels_AUS_attribute_data.lc08_rfcropp;
CAMELS_AUS_data.lc09_rfpastu = camels_AUS_attribute_data.lc09_rfpastu;
CAMELS_AUS_data.lc10_rfsugar = camels_AUS_attribute_data.lc10_rfsugar;
CAMELS_AUS_data.lc11_wetlands = camels_AUS_attribute_data.lc11_wetlands;
CAMELS_AUS_data.lc14_tussclo = camels_AUS_attribute_data.lc14_tussclo;
CAMELS_AUS_data.lc15_alpineg = camels_AUS_attribute_data.lc15_alpineg;
CAMELS_AUS_data.lc16_openhum = camels_AUS_attribute_data.lc16_openhum;
CAMELS_AUS_data.lc18_opentus = camels_AUS_attribute_data.lc18_opentus;
CAMELS_AUS_data.lc19_shrbsca = camels_AUS_attribute_data.lc19_shrbsca;
CAMELS_AUS_data.lc24_shrbden = camels_AUS_attribute_data.lc24_shrbden;
CAMELS_AUS_data.lc25_shrbope = camels_AUS_attribute_data.lc25_shrbope;
CAMELS_AUS_data.lc31_forclos = camels_AUS_attribute_data.lc31_forclos;
CAMELS_AUS_data.lc32_foropen = camels_AUS_attribute_data.lc32_foropen;
CAMELS_AUS_data.lc33_woodope = camels_AUS_attribute_data.lc33_woodope;
CAMELS_AUS_data.lc34_woodspa = camels_AUS_attribute_data.lc34_woodspa;
CAMELS_AUS_data.lc35_urbanar = camels_AUS_attribute_data.lc35_urbanar;
CAMELS_AUS_data.prop_forested = camels_AUS_attribute_data.prop_forested;
CAMELS_AUS_data.nvis_grasses_n = camels_AUS_attribute_data.nvis_grasses_n;
CAMELS_AUS_data.nvis_grasses_e = camels_AUS_attribute_data.nvis_grasses_e;
CAMELS_AUS_data.nvis_forests_n = camels_AUS_attribute_data.nvis_forests_n;
CAMELS_AUS_data.nvis_forests_e = camels_AUS_attribute_data.nvis_forests_e;
CAMELS_AUS_data.nvis_shrubs_n = camels_AUS_attribute_data.nvis_shrubs_n;
CAMELS_AUS_data.nvis_shrubs_e = camels_AUS_attribute_data.nvis_shrubs_e;
CAMELS_AUS_data.nvis_woodlands_n = camels_AUS_attribute_data.nvis_woodlands_n;
CAMELS_AUS_data.nvis_woodlands_e = camels_AUS_attribute_data.nvis_woodlands_e;
CAMELS_AUS_data.nvis_bare_n = camels_AUS_attribute_data.nvis_bare_n;
CAMELS_AUS_data.nvis_bare_e = camels_AUS_attribute_data.nvis_bare_e;
CAMELS_AUS_data.nvis_nodata_n = camels_AUS_attribute_data.nvis_nodata_n;
CAMELS_AUS_data.nvis_nodata_e = camels_AUS_attribute_data.nvis_nodata_e;

% anthropogenic influences						
CAMELS_AUS_data.distupdamw = camels_AUS_attribute_data.distupdamw;
CAMELS_AUS_data.impound_fac = camels_AUS_attribute_data.impound_fac;
CAMELS_AUS_data.flow_div_fac = camels_AUS_attribute_data.flow_div_fac;
CAMELS_AUS_data.leveebank_fac = camels_AUS_attribute_data.leveebank_fac;
CAMELS_AUS_data.infrastruc_fac = camels_AUS_attribute_data.infrastruc_fac;
CAMELS_AUS_data.settlement_fac = camels_AUS_attribute_data.settlement_fac;
CAMELS_AUS_data.extract_ind_fac = camels_AUS_attribute_data.extract_ind_fac;
CAMELS_AUS_data.landuse_fac = camels_AUS_attribute_data.landuse_fac;
CAMELS_AUS_data.catchment_di = camels_AUS_attribute_data.catchment_di;
CAMELS_AUS_data.flow_regime_di = camels_AUS_attribute_data.flow_regime_di;
CAMELS_AUS_data.river_di = camels_AUS_attribute_data.river_di;

% other 
CAMELS_AUS_data.pop_mean = camels_AUS_attribute_data.pop_mean;
CAMELS_AUS_data.pop_max = camels_AUS_attribute_data.pop_max;
CAMELS_AUS_data.pop_gt_1 = camels_AUS_attribute_data.pop_gt_1;
CAMELS_AUS_data.pop_gt_10 = camels_AUS_attribute_data.pop_gt_10;
CAMELS_AUS_data.erosivity = camels_AUS_attribute_data.erosivity;
CAMELS_AUS_data.anngro_mega = camels_AUS_attribute_data.anngro_mega;
CAMELS_AUS_data.anngro_meso = camels_AUS_attribute_data.anngro_meso;
CAMELS_AUS_data.anngro_micro = camels_AUS_attribute_data.anngro_micro;
CAMELS_AUS_data.gromega_seas = camels_AUS_attribute_data.gromega_seas;
CAMELS_AUS_data.gromeso_seas = camels_AUS_attribute_data.gromeso_seas;
CAMELS_AUS_data.gromicro_seas = camels_AUS_attribute_data.gromicro_seas;
CAMELS_AUS_data.npp_ann = camels_AUS_attribute_data.npp_ann;
CAMELS_AUS_data.npp_1 = camels_AUS_attribute_data.npp_1;
CAMELS_AUS_data.npp_2 = camels_AUS_attribute_data.npp_2;
CAMELS_AUS_data.npp_3 = camels_AUS_attribute_data.npp_3;
CAMELS_AUS_data.npp_4 = camels_AUS_attribute_data.npp_4;
CAMELS_AUS_data.npp_5 = camels_AUS_attribute_data.npp_5;
CAMELS_AUS_data.npp_5 = camels_AUS_attribute_data.npp_5;
CAMELS_AUS_data.npp_6 = camels_AUS_attribute_data.npp_6;
CAMELS_AUS_data.npp_7 = camels_AUS_attribute_data.npp_7;
CAMELS_AUS_data.nnp_8 = camels_AUS_attribute_data.npp_8;
CAMELS_AUS_data.npp_9 = camels_AUS_attribute_data.npp_9;
CAMELS_AUS_data.npp_10 = camels_AUS_attribute_data.npp_10;
CAMELS_AUS_data.npp_11 = camels_AUS_attribute_data.npp_11;
CAMELS_AUS_data.npp_12 = camels_AUS_attribute_data.npp_12;

%% Load hydro-meteorological time series
% To extract the time series, we loop over all catchments. We also
% calculate the completeness of the flow records.
flow_perc_complete = NaN(length(CAMELS_AUS_data.station_id),1);
P = cell(length(CAMELS_AUS_data.station_id),1); % precipitation
PET = cell(length(CAMELS_AUS_data.station_id),1); % potential evapotranspiration
Q = cell(length(CAMELS_AUS_data.station_id),1); % streamflow
T = cell(length(CAMELS_AUS_data.station_id),1); % temperature

Q_table = readtable(...
    strcat(path,'03_streamflow/streamflow_mmd'),'ReadVariableNames',true);
Q_date = datenum(table2array(Q_table(:,1:3)));

P_table = readtable(...
    strcat(path,'05_hydrometeorology/01_precipitation_timeseries/precipitation_AWAP'),...
    'ReadVariableNames',true);
P_date = datenum(table2array(P_table(:,1:3)));

PET_table = readtable(...
    strcat(path,'05_hydrometeorology/02_EvaporativeDemand_timeseries/et_morton_wet_SILO'),...
    'ReadVariableNames',true);
PET_date = datenum(table2array(PET_table(:,1:3)));

Tmax_table = readtable(...
    strcat(path,'05_hydrometeorology/03_Other/AWAP/tmax_AWAP'),'ReadVariableNames',true);
Tmin_table = readtable(...
    strcat(path,'05_hydrometeorology/03_Other/AWAP/tmin_AWAP'),'ReadVariableNames',true);
T_date = datenum(table2array(Tmax_table(:,1:3)));
% Tmin_date = datenum(table2array(Tmin_table(:,1:3)));

fprintf('Loading catchment data (AUS)...\n')
for i = 1:size(Q_table,2)-3 % loop over all catchments
    
    if mod(i,100) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,length(CAMELS_AUS_data.station_id))
    end
    
    [P{i}, PET{i}, Q{i}, T{i}] = loadCatchmentCAMELS_AUS(...
        i+3,Q_table,Q_date,P_table,P_date,PET_table,PET_date,...
        Tmax_table,Tmin_table,T_date);
    flow_perc_complete(i) = 100*(1-sum(isnan(Q{i}(:,2)))./length(Q{i}(:,2)));
    
end

% add hydro-meteorological time series to struct file
CAMELS_AUS_data.flow_perc_complete = flow_perc_complete;
CAMELS_AUS_data.P = P;
CAMELS_AUS_data.PET = PET;
CAMELS_AUS_data.Q = Q;
CAMELS_AUS_data.T = T;

% save the struct file
if save_struct
    save('CAMELS_Matlab/Data/CAMELS_AUS_data.mat','-struct','CAMELS_AUS_data')
end

end
