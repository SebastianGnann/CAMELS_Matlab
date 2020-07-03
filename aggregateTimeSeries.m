function [X_annual, X_monthly, year_list] = ...
    aggregateTimeSeries(X, start_water_year)
%aggregateTimeSeries Calculate annual and monthly data from time series.
%
%   INPUT
%   X: time series (datenum, values)
%   start_water_year: first month of water year (default: 1, i.e. January)
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
    start_water_year = 1;
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
    X_water_year = ...
        [X_temp(years==year & months>=start_water_year); ...
        X_temp(years==year+1 & months<start_water_year)];
    X_annual(y) = mean(X_water_year);
    % check again
    for m = start_water_year:12
        X_monthly(y,m-start_water_year+1) = mean(X_temp(years==year & months==m));
    end
    for m = 1:start_water_year-1
        X_monthly(y,m+13-start_water_year) = mean(X_temp(years==year+1 & months==m));
    end
end

end