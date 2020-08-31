function [P, PET, Q, T] = loadCatchmentCAMELS(ID,path)
%loadCatchmentCAMELS Loads catchment data (P, PET, Q, T)
%
%   INPUT
%   ID: catchment ID
%   path: file path
%
%   OUTPUT
%   P: precipitation [mm/d]
%   PET: potential evapotranspiration (adjusted using standard coefficient
%   of 1.26) [mm/d]
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

% path = 'C:\Users\sg16200\Local Documents\Camels\basin_timeseries_v1p2_modelOutput_daymet\model_output_daymet\model_output\flow_timeseries\daymet\_all\';
foundID = false; % search through all folders until catchment is found
i_str_list = ["01","02","03","04","05","06","07","08","09","10","11",...
    "12","13","14","15","16","17","18"];
i = 0;
while foundID == false
    i = i+1;
    i_str = i_str_list(i);
    try
        file_ID_model = strcat(path,i_str,'\',num2str(ID,'%08d'),'_05_model_output.txt');
        txt_data=fileread(file_ID_model);
        file_ID_parameters = strcat(path,i_str,'\',num2str(ID,'%08d'),'_05_model_parameters.txt');
        txt_para=fileread(file_ID_parameters);
        foundID = true;
    catch
        disp('')
        foundID = false;
    end
end

% NOTE: NOW version 05 of Newman dataset
% YR MNTH DY HR SWE PRCP RAIM TAIR PET ET MOD_RUN OBS_RUN

data_model_cell = textscan(txt_data,...
    '%f %f %f %f %f %f %f %f %f %f %f %f', ...
    'Delimiter', '\t', 'HeaderLines', 1);

data_parameter_cell = textscan(txt_para,...
    '%s %f', ...
    'Delimiter', '\t', 'HeaderLines', 0);
PET_coefficient = data_parameter_cell{2}(41);
% disp(PET_coefficient)

Y = data_model_cell{1};
M = data_model_cell{2};
D = data_model_cell{3};
date = datenum(Y,M,D);

Q_temp = data_model_cell{12};
Q_temp(Q_temp==-999) = NaN;
P_temp = data_model_cell{6};
PET_temp = data_model_cell{9}; 
T_temp = data_model_cell{8};

Q = [date Q_temp];
P = [date P_temp];
% PET = [date PET_temp]; % calibrated PET
PET = [date (1.26/PET_coefficient).*PET_temp]; % adjusted PET
T = [date T_temp];

end

