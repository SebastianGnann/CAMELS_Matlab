function [P, PET, Q, T] = loadCatchmentCAMELSGB(ID,path)
%loadCatchmentCAMELSGB Loads catchment data (P, PET, Q, T) and creates 
%   name-strings (for CAMELS-GB format)
%
%   INPUT
%   ID: string with station ID (or catchment ID)
%   path: file path
%
%   OUTPUT
%   P: precipitation [mm/d]
%   PET: potential evapotranspiration (without interception correction) [mm/d]
%   Q: streamflow [mm/d]
%   T: T [°C]
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

% check input parameters
if nargin < 2
    error('Not enough input arguments.')
end

% path = "C:\Users\sg16200\Local Documents\CAMELS_GB\CAMELS_GB_Timeseries\";

file_ID = strcat(path,'CAMELS_GB_hydromet_timeseries_',num2str(ID),'_19701001-20150930.csv');

% date	precipitation	pet	temperature	discharge_spec	discharge_vol	peti	humidity	shortwave_rad	longwave_rad	windspeed
[data,data_str] = xlsread(file_ID);

formatIn = 'dd/mm/yyyy';
date = datenum(data_str(2:end,1),formatIn);

Q_temp = data(:,4);
P_temp = data(:,1);
PET_temp = data(:,2);
% PET_temp = data(:,6); % with interception correction
T_temp = data(:,3);

Q = [date Q_temp];
P = [date P_temp];
% PET = [date PET_temp];
PET = [date PET_temp];
T = [date T_temp];

%{ 
txt=fileread(file_ID);
% YYYY	MM	DD	RAIN	PET	TEMPERATURE	DISCHARGE_MM	DISCHARGE_M3S	PETI	HUMIDITY	SWR	LWR	WINDSPEED
data_cell = textscan(txt,...
    '%f %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'Delimiter', '\t', 'HeaderLines', 1);

% Y = data_cell{1};
% M = data_cell{2};
% D = data_cell{3};
% date = datenum(Y,M,D);

% Q_temp = data_cell{7};
% % Q_temp = data_cell{8};
% Q_temp(Q_temp==-999) = NaN;
% P_temp = data_cell{4};
% PET_temp = data_cell{5};
% % PET_temp = data_cell{9}; % with interception correction
% T_temp = data_cell{6};
% 
% Q = [date Q_temp];
% P = [date P_temp];
% % PET = [date PET_temp];
% PET = [date PET_temp];
% T = [date T_temp];
%}

end

