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
t_string = datetime(t,'ConvertFrom','datenum');
[years, months, days] = ymd(t_string);

if start_water_year == 1
    % calendar year
    year_start = min(years);
    year_end = max(years);
else
    % water year always corresponds to the last day of the water year,
    % e.g. water year starting from 1 October 1999 is water year 2000
    year_start = min(years)+1;
    year_end = max(years);
end
year_list = [year_start:year_end]';

if months(1) ~= start_water_year && days(1)~=1
    warning('Time series and chosen water year do not match. Incomplete years possible.')
end
if months(end) ~= start_water_year-1 && days(end)<28
    warning('Time series and chosen water year do not match. Incomplete years possible.')
end

X_annual = NaN(year_end-year_start,1);
X_monthly = NaN(year_end-year_start,12);

% extract years, months, etc.
for y = 1:length(year_list)
    year = year_list(y);
    X_water_year = ...
        [X(years==year-1 & months>=start_water_year); ...
        X(years==year & months<start_water_year)];
    X_annual(y) = sum(X_water_year,'omitnan');
    % check again
    for m = start_water_year:12
        X_monthly(y,m-start_water_year+1) = sum(X(years==year-1 & months==m),'omitnan');
    end
    for m = 1:start_water_year-1
        X_monthly(y,m+13-start_water_year) = sum(X(years==year & months==m),'omitnan');
    end
end

end