function [P, PET, Q, T] = loadCatchmentCAMELSBR(ID,path)
%loadCatchmentCAMELSBR Loads catchment data (P, PET, Q, T) and creates
%   name-strings (for CAMELS BR format)
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

% streamflow
file_ID_model = strcat(path,'3_CAMELS_BR_streamflow_mm','\',num2str(ID,'%08d'),'_streamflow_mm.txt');
txt_data=fileread(file_ID_model);
data_model_cell = textscan(txt_data,'%f %f %f %f %f %f','Delimiter',' ', 'HeaderLines', 1);
Y = data_model_cell{1};
M = data_model_cell{2};
D = data_model_cell{3};
date = datenum(Y,M,D);
Q_temp = data_model_cell{4};
Q = [date Q_temp];

% precipitation CHIRPS
file_ID_model = strcat(path,'4_CAMELS_BR_precipitation_chirps','\',num2str(ID,'%08d'),'_precipitation_chirps.txt');
txt_data=fileread(file_ID_model);
data_model_cell = textscan(txt_data,'%f %f %f %f','Delimiter',' ', 'HeaderLines', 1);
Y = data_model_cell{1};
M = data_model_cell{2};
D = data_model_cell{3};
date = datenum(Y,M,D);
P_temp = data_model_cell{4};
P = [date P_temp];

% pot. evapotranspiration GLEAM
file_ID_model = strcat(path,'7_CAMELS_BR_evapotransp_gleam','\',num2str(ID,'%08d'),'_evapotransp_gleam.txt');
txt_data=fileread(file_ID_model);
data_model_cell = textscan(txt_data,'%f %f %f %f','Delimiter',' ', 'HeaderLines', 1);
Y = data_model_cell{1};
M = data_model_cell{2};
D = data_model_cell{3};
date = datenum(Y,M,D);
PET_temp = data_model_cell{4};
PET = [date PET_temp];

% temperature CPC
file_ID_model = strcat(path,'11_CAMELS_BR_temperature_mean_cpc','\',num2str(ID,'%08d'),'_temperature_mean.txt');
txt_data=fileread(file_ID_model);
data_model_cell = textscan(txt_data,'%f %f %f %f','Delimiter',' ', 'HeaderLines', 1);
Y = data_model_cell{1};
M = data_model_cell{2};
D = data_model_cell{3};
date = datenum(Y,M,D);
T_temp = data_model_cell{4};
T = [date T_temp];

end

