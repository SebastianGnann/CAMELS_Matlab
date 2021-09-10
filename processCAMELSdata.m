%% processCAMELSdata
% 
%   This script shows how to process the CAMELS struct files created with
%   saveCAMELSdata.m; specifically, it shows how to:
%   - load the struct file,
%   - use the data to calculate some metrics (using the TOSSH toolbox), 
%   - plot maps of those metrics (requires the BrewerMap toolbox).
%
%   Note that this script requires sufficient RAM.
%
%   Copyright (C) 2021
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

%% load useful packages 

if (exist('BrewerMap') == 7)
    addpath(genpath('BrewerMap'));
else
    error('BrewerMap toolbox needed. Can be downloaded from https://github.com/DrosteEffect/BrewerMap and should be in a folder named BrewerMap in the same directory.')
end

if (exist('TOSSH') == 7)
    addpath(genpath('TOSSH'));
else
    error('TOSSH toolbox needed. Can be downloaded from https://github.com/TOSSHtoolbox and should be in a folder named TOSSH in the same directory.')
end

%% add paths

% working directory (important so that functions herein are called)
mydir = 'CAMELS_Matlab';
addpath(genpath(mydir));

% figure path
fig_path = 'CAMELS_Matlab/Figures/';
results_path = 'CAMELS_Matlab/Results/';

%% load catchment data

CAMELS_US_data = load('CAMELS_Matlab/Data/CAMELS_US_data.mat');
% CAMELS_GB_data = load('CAMELS_Matlab/Data/CAMELS_GB_data.mat');
% CAMELS_CL_data = load('CAMELS_Matlab/Data/CAMELS_CL_data.mat');
% CAMELS_BR_data = load('CAMELS_Matlab/Data/CAMELS_BR_data.mat');
% CAMELS_AUS_data = load('CAMELS_Matlab/Data/CAMELS_AUS_data.mat');

%% Calculate signatures using TOSSH
% To use TOSSH calculation functions, we need to create cell arrays
% containing the time series. We use cell arrays since some time series
% might have different lengths. While the length of each row in the cell
% array can vary, the cell arrays containing the t, Q, P, and PET data need
% to have exactly the same dimensions. We first initialise the cell arrays.
n_CAMELS_US = length(CAMELS_US_data.gauge_id);
t_mat_US = cell(n_CAMELS_US,1);
Q_mat_US = cell(n_CAMELS_US,1);
P_mat_US = cell(n_CAMELS_US,1);
PET_mat_US = cell(n_CAMELS_US,1);

% We then loop over all catchments and extract the time series for each
% catchment.
fprintf('Creating data matrix...\n')
for i = 1:n_CAMELS_US
    
    if mod(i,100) == 0 % check progress
        fprintf('%.0f/%.0f\n',i,n_CAMELS_US)
    end
    
    t = datetime(CAMELS_US_data.Q{i}(:,1),'ConvertFrom','datenum');
    Q = CAMELS_US_data.Q{i}(:,2);    
    P = CAMELS_US_data.P{i}(:,2);
    PET = CAMELS_US_data.PET{i}(:,2);
    
    % get subperiod, here from 1 Oct 1989 to 30 Sep 2009
    indices = 1:length(t); 
    start_ind = indices(t==datetime(1989,10,1));
    % in case time series starts after 1 Oct 1989
    if isempty(start_ind); start_ind = 1; end 
    end_ind = indices(t==datetime(2009,9,30));    
    t = t(start_ind:end_ind);
    Q = Q(start_ind:end_ind);
    P = P(start_ind:end_ind);
    PET = PET(start_ind:end_ind);
    % calculate completeness during sub-period
    flow_perc_complete = 100*(1-sum(isnan(Q))./length(Q));
    CAMELS_data.flow_perc_complete(i) = flow_perc_complete;
    
    t_mat_US{i} = t;
    Q_mat_US{i} = Q;
    P_mat_US{i} = P;
    PET_mat_US{i} = PET;
    
end

% We can now use the calculation functions to calculate various signatures.
% McMillan groundwater set - note that this make take a few minutes
CAMELS_US_signatures_Groundwater = ...
    calc_McMillan_Groundwater(Q_mat_US, t_mat_US, P_mat_US, PET_mat_US);

% We can save the results as mat files which can be easily loaded. 
save(strcat(results_path,'CAMELS_US_signatures_Groundwater.mat'),...
    '-struct','CAMELS_US_signatures_Groundwater')

% Alternatively, we can save the results as txt file.
% writetable(struct2table(CAMELS_US_signatures_Groundwater),...
%     strcat(results_path,'CAMELS_US_signatures_Groundwater.txt'))

%% correlation matrix
Groundwater_matrix = [...
    CAMELS_US_signatures_Groundwater.AverageStorage,...
    CAMELS_US_signatures_Groundwater.BaseflowRecessionK,...
    CAMELS_US_signatures_Groundwater.BFI,...
    CAMELS_US_signatures_Groundwater.EventRR,...
    CAMELS_US_signatures_Groundwater.EventRR_TotalRR_ratio,...
    CAMELS_US_signatures_Groundwater.MRC_num_segments,...
    CAMELS_US_signatures_Groundwater.Recession_a_Seasonality,...
    CAMELS_US_signatures_Groundwater.RecessionParameters(:,1),...
    CAMELS_US_signatures_Groundwater.RecessionParameters(:,2),...
    CAMELS_US_signatures_Groundwater.RR_Seasonality,...
    CAMELS_US_signatures_Groundwater.Spearmans_rho,...
    CAMELS_US_signatures_Groundwater.StorageFraction(:,1),...
    CAMELS_US_signatures_Groundwater.TotalRR,...
    CAMELS_US_signatures_Groundwater.VariabilityIndex,...
    ];

signature_names = {...
    'AverageStorage',...
    'BaseflowRecessionK',...
    'BFI',...
    'EventRR',...
    'EventRR_TotalRR_ratio',...
    'MRC_num_segments',...
    'Recession_a_Seasonality',...
    'RecessionParameters_a',...
    'RecessionParameters_b',...
    'RR_Seasonality',...
    'Spearmans_rho',...
    'StorageFraction',...
    'TotalRR',...
    'VariabilityIndex',...
    };

rho = corr(Groundwater_matrix,'rows','complete');
figure('pos',[100 100 380 300]);
imagesc(rho); % axis equal
set(gca, 'XTick', 1:length(rho)); % center x-axis ticks on bins
set(gca, 'YTick', 1:length(rho)); % center y-axis ticks on bins
% set(gca, 'XTickLabel', signature_names); % set x-axis labels
% set(gca, 'YTickLabel', signature_names); % set y-axis labels
colormap('parula'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1])

%% plot results
% We can plot and save maps of some attribute using the plotMap functions.

% Groundwater set
plotMap_US(CAMELS_US_data.gauge_lat,CAMELS_US_data.gauge_lon,...
    CAMELS_US_signatures_Groundwater.AverageStorage,...
    'attribute_name','Av. storage [mm]',...
    'ID',CAMELS_US_data.gauge_id,...
    'colour_scheme','RdBu','flip_colour_scheme',true,...
    'c_limits',[0 500],'c_log_scale',false,...
    'figure_name','AverageStorage','save_figure',false,'figure_path',fig_path)

plotMap_US(CAMELS_US_data.gauge_lat,CAMELS_US_data.gauge_lon,...
    1./CAMELS_US_signatures_Groundwater.BaseflowRecessionK,...
    'attribute_name','Rec. K [days]',...
    'ID',CAMELS_US_data.gauge_id,...
    'colour_scheme','YlGnBu','flip_colour_scheme',false,...
    'c_limits',[0 30],'c_log_scale',false,...
    'figure_name','RecessionK','save_figure',false,'figure_path',fig_path)

plotMap_US(CAMELS_US_data.gauge_lat,CAMELS_US_data.gauge_lon,...
    CAMELS_US_signatures_Groundwater.BFI,...
    'attribute_name','BFI [-]',...
    'ID',CAMELS_US_data.gauge_id,...
    'colour_scheme','Spectral','flip_colour_scheme',false,...
    'c_limits',[0 1],'c_log_scale',false,...
    'figure_name','BFI','save_figure',false,'figure_path',fig_path)

