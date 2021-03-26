%% useCAMELSdata
% 
%   This script shows how to use the CAMELS struct files created with
%   saveCAMELSdata. Specifically, it shows how to:
%   - load the files 
%   - plot the data (scatter plots, maps)
%   - use the time series to calculate some metrics. 
%
%   Note that this script requires sufficient RAM.
%
%   Copyright (C) 2021
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

%% add paths

% working directory
mydir = 'CAMELS_Matlab';
addpath(genpath(mydir));

% figure path
fig_path = 'CAMELS_Matlab/Figures';

%% load useful packages 

if exist('BrewerMap') == 7
    addpath(genpath('BrewerMap'));
else
    error('BrewerMap toolbox needed. Can be downloaded from https://github.com/DrosteEffect/BrewerMap and should be in a folder named BrewerMap in the same directory.')
end

%% load catchment data

CAMELS_US_data = load('CAMELS_Matlab/Data/CAMELS_US_data.mat');
CAMELS_GB_data = load('CAMELS_Matlab/Data/CAMELS_GB_data.mat');
CAMELS_CL_data = load('CAMELS_Matlab/Data/CAMELS_CL_data.mat');
CAMELS_BR_data = load('CAMELS_Matlab/Data/CAMELS_BR_data.mat');
CAMELS_AUS_data = load('CAMELS_Matlab/Data/CAMELS_AUS_data.mat');

%% prepare data

subsetUS = true(length(CAMELS_US_data.gauge_id),1);
subsetUS(CAMELS_US_data.frac_snow>0.2) = false;

subsetGB = true(length(CAMELS_GB_data.gauge_id),1);
subsetGB(CAMELS_GB_data.frac_snow>0.2) = false;
% subsetGB(~CAMELS_GB_data.isBenchmark) = false;

subsetCL = true(length(CAMELS_CL_data.gauge_id),1);
% subsetCL(CAMELS_CL_data.interv_degree>1) = false;
subsetCL(CAMELS_CL_data.frac_snow_mswep>0.2) = false;
subsetCL(CAMELS_CL_data.lc_glacier>1) = false;

subsetBR = true(length(CAMELS_BR_data.gauge_id),1);
subsetBR(CAMELS_BR_data.frac_snow>0.2) = false;
% subsetBR(CAMELS_BR_data.degree_of_regulation>1) = false;

subsetAUS = true(length(CAMELS_AUS_data.station_id),1);
subsetAUS(CAMELS_AUS_data.frac_snow>0.2) = false;


%% to do

% - maps
% - scatter plots


%% scatter plots

% We can have a look at the different brewermap color scales by running:
% brewermap('plot')

% scatter plots
figure; hold on
scatter(CAMELS_US_data.aridity(subsetUS),CAMELS_US_data.runoff_ratio(subsetUS),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
scatter(CAMELS_GB_data.aridity(subsetGB),CAMELS_GB_data.runoff_ratio(subsetGB),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
scatter(CAMELS_CL_data.aridity_mswep(subsetCL),CAMELS_CL_data.runoff_ratio_mswep(subsetCL),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
scatter(CAMELS_BR_data.aridity(subsetBR),CAMELS_BR_data.runoff_ratio(subsetBR),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
scatter(CAMELS_AUS_data.aridity(subsetAUS),CAMELS_AUS_data.runoff_ratio(subsetAUS),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
legend('US','GB','CL','BR','AU')
xlabel('PET/P'); ylabel('Q/P')
xlim([.1 10]); ylim([0 1])
set(gca,'xscale','log')

figure; hold on
scatter(CAMELS_US_data.aridity(subsetUS),CAMELS_US_data.baseflow_index(subsetUS).*CAMELS_US_data.runoff_ratio(subsetUS),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
scatter(CAMELS_GB_data.aridity(subsetGB),CAMELS_GB_data.baseflow_index(subsetGB).*CAMELS_GB_data.runoff_ratio(subsetGB),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
scatter(CAMELS_CL_data.aridity_mswep(subsetCL),CAMELS_CL_data.baseflow_index(subsetCL).*CAMELS_CL_data.runoff_ratio_mswep(subsetCL),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
scatter(CAMELS_BR_data.aridity(subsetBR),CAMELS_BR_data.baseflow_index(subsetBR).*CAMELS_BR_data.runoff_ratio(subsetBR),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
scatter(CAMELS_AUS_data.aridity(subsetAUS),CAMELS_AUS_data.baseflow_index(subsetAUS).*CAMELS_AUS_data.runoff_ratio(subsetAUS),25,'filled','markerfacealpha',.25,'markerfacecolor','k')
legend('US','GB','CL','BR','AU')
xlabel('PET/P'); ylabel('Qb/P')
xlim([.1 10]); ylim([0 .7])
set(gca,'xscale','log')

%% map plots

% We can plot a map of some attribute using the plotMap functions. 
plotMapUS(CAMELS_US_data.gauge_lat,CAMELS_US_data.gauge_lon,...
    CAMELS_US_data.aridity,'save_figure',false)

% Various options can be specified and the figure can be saved.
plotMapUS(CAMELS_US_data.gauge_lat,CAMELS_US_data.gauge_lon,...
    CAMELS_US_data.aridity,'attribute_name','PET/P [-]',...
    'ID',CAMELS_US_data.gauge_id,...
    'colour_scheme','RdBu','flip_colour_scheme',true,...
    'c_limits',[0.1 10],'c_log_scale',true,...
    'c_lower_limit_open',false,'c_upper_limit_open',false,...
    'figure_title','','figure_name','Aridity',...
    'save_figure',false,'figure_path',fig_path,'figure_type','-dpng')

% We can plot maps for various countries.
% to do: adjust colorbars

plotMapUS(CAMELS_US_data.gauge_lat,CAMELS_US_data.gauge_lon,...
    CAMELS_US_data.aridity,'attribute_name','PET/P [-]',...
    'ID',CAMELS_US_data.gauge_id,...
    'colour_scheme','RdBu','flip_colour_scheme',true,...
    'c_limits',[0.1 10],'c_log_scale',true)

plotMapGB(CAMELS_GB_data.gauge_lat,CAMELS_GB_data.gauge_lon,...
    CAMELS_GB_data.aridity,'attribute_name','PET/P [-]',...
    'ID',CAMELS_GB_data.gauge_id,...
    'colour_scheme','RdBu','flip_colour_scheme',true,...
    'c_limits',[0.1 10],'c_log_scale',true)

plotMapCL(CAMELS_CL_data.gauge_lat,CAMELS_CL_data.gauge_lon,...
    CAMELS_CL_data.aridity_mswep,'attribute_name','PET/P [-]',...
    'ID',CAMELS_CL_data.gauge_id,...
    'colour_scheme','RdBu','flip_colour_scheme',true,...
    'c_limits',[0.1 10],'c_log_scale',true)

plotMapBR(CAMELS_BR_data.gauge_lat,CAMELS_BR_data.gauge_lon,...
    CAMELS_BR_data.aridity,'attribute_name','PET/P [-]',...
    'ID',CAMELS_BR_data.gauge_id,...
    'colour_scheme','RdBu','flip_colour_scheme',true,...
    'c_limits',[0.1 10],'c_log_scale',true)

plotMapAUS(CAMELS_AUS_data.lat_centroid,CAMELS_AUS_data.long_centroid,...
    CAMELS_AUS_data.aridity,'attribute_name','PET/P [-]',...
    'ID',CAMELS_AUS_data.gauge_id,...
    'colour_scheme','RdBu','flip_colour_scheme',true,...
    'c_limits',[0.1 10],'c_log_scale',true)

%% new signatures

