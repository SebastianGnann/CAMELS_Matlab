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
%   PET: potential evapotranspiration [mm/d]
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

file_ID_model = strcat(path,'CAMELS_GB_hydromet_timeseries_',num2str(ID),'.txt');
txt=fileread(file_ID_model);

% YYYY	MM	DD	RAIN	PET	TEMPERATURE	DISCHARGE_MM	DISCHARGE_M3S	PETI	HUMIDITY	SWR	LWR	WINDSPEED

data_model_cell = textscan(txt,...
    '%f %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'Delimiter', '\t', 'HeaderLines', 1);

Y = data_model_cell{1};
M = data_model_cell{2};
D = data_model_cell{3};
date = datenum(Y,M,D);

Q_temp = data_model_cell{7};
% Q_temp = data_model_cell{8};
Q_temp(Q_temp==-999) = NaN;
P_temp = data_model_cell{4};
% PET_temp = data_model_cell{4};
PET_temp = data_model_cell{9}; % with interception correction
T_temp = data_model_cell{5};

Q = [date Q_temp];
P = [date P_temp];
PET = [date PET_temp];
T = [date T_temp];

end

