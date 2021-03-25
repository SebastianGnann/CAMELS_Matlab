function [P, PET, Q, T] = loadCatchmentCAMELS_CL(...
    index,Q_table,Q_date,P_table,P_date,PET_table,PET_date,T_table,T_date)
%loadCatchmentCAMELS_CL Loads catchment data (P, PET, Q, T) for CAMELS CL 
%   format.
%
%   INPUT
%   index: index of catchment (1 is first catchment in list)
%   X_table: tables with data
%   X_date: vectors with dates
%
%   OUTPUT
%   P: precipitation [mm/d]
%   PET: potential evapotranspiration [mm/d]
%   Q: streamflow [mm/d]
%   T: T [°C]
%
%   Copyright (C) 2021
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 2
    error('Not enough input arguments.')
end

Q_temp = table2array(Q_table(2:end,index+1));
P_temp = table2array(P_table(2:end,index+1));
PET_temp = table2array(PET_table(2:end,index+1));
T_temp = table2array(T_table(2:end,index+1));

% get earliest/latest common date to get consistent time series
date_min = max([min(Q_date);min(P_date);min(PET_date);min(T_date)]);
date_max = min([max(Q_date);max(P_date);max(PET_date);max(T_date)]);

ind = 1:length(Q_date);
Q_range = [ind(Q_date==date_min):ind(Q_date==date_max)];
ind = 1:length(P_date);
P_range = [ind(P_date==date_min):ind(P_date==date_max)];
ind = 1:length(PET_date);
PET_range = [ind(PET_date==date_min):ind(PET_date==date_max)];
ind = 1:length(T_date);
T_range = [ind(T_date==date_min):ind(T_date==date_max)];

Q = [Q_date(Q_range) Q_temp(Q_range)];
P = [P_date(P_range) P_temp(P_range)];
PET = [PET_date(PET_range) PET_temp(PET_range)];
T = [T_date(T_range) T_temp(T_range)];

end

