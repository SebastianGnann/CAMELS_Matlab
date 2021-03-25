function [X] = getSubPeriod(date_start,date_end,X)
%getSubPeriod Extract subperiod from time series, e.g. from 
%   1 Oct 1989 to 30 Sept 1999.
%
%   INPUT
%   date_start: start date as Matlab datenum
%   date_end: end date as Matlab datenum
% 	X: time series (datenum,value)
%
%   OUTPUT
%   trimmed time series (datenum,value)
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

if isempty(X)
    % ignore empty cells
else
    date_temp = X(:,1);
    index_start = find(date_temp == date_start);
    index_end = find(date_temp == date_end);
    X = [X(index_start:index_end,1) X(index_start:index_end,2)];
end

end



