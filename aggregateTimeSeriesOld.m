function [X_annual, X_monthly, year_list] = ...
    aggregateTimeSeriesOld(X, useWaterYear) 
%aggregateTimeSeriesOld Calculate annual and monthly data from time series.
%
%   INPUT
%   X: time series (datenum, values)
%   useWaterYear: use water year instead (starting from 1st Oct) (logical)
%
%   OUTPUT
%   X_annual: annual values
%   X_monthly: monthly values
%   year_list: corresponding years
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

if nargin < 2
    useWaterYear = false;
end

% get years and months
t_string = datetime(X(:,1),'ConvertFrom','datenum');
[years, months, ~] = ymd(t_string);
year_start = min(years);
year_end = max(years);

year_list = [year_start:year_end]';
% water_year = NaN(year_end-year_start,1);

X_temp = X(:,2);

X_annual = NaN(year_end-year_start,1);
X_monthly = NaN(year_end-year_start,12);

% extract years, months, etc.
for y = 1:length(year_list)    
    year = year_list(y);    
    if useWaterYear        
        X_water_year = ...
            [X_temp(years==year & months>=10); ...
            X_temp(years==year+1 & months<10)];
        X_annual(y) = mean(X_water_year);    
        % check again
        for m = 10:12
            X_monthly(y,m-9) = mean(X_temp(years==year & months==m));
        end             
        for m = 1:9
            X_monthly(y,m+3) = mean(X_temp(years==year+1 & months==m));
        end 
    else
        X_year = X_temp(years==year);
        X_annual(y) = mean(X_year);              
        for m = 1:12
            X_monthly(y,m) = mean(X_temp(years==year & months==m));            
        end        
    end
end

end