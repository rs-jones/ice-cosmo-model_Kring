%% CREATE TIME-SERIES FILE FOR COMPLEX HISTORY MODEL
%
% Loads time-series data from an Excel spreadsheet (.xlsx), and creates a
% new file with necessary information.
%
% Input data must be 2 columns (time (years before present), values).
%
%
%% Enter Information

clear; % Start fresh

% SET (shown with examples)
data_name = 'SST_SO_stack';                 % Name used for .xlsx
ts_name = 'Southern Ocean sea surface temperature'; % Name of time series
ts_short_name = 'SST_{SO}';                 % Short name of time series
ts_unit = '°C';                             % Unit of time-series data


%% Run

% Get data
[in_data,~,~] = xlsread(strcat(data_name,'.xlsx'));

% Check data
if (numel(in_data(1,:)) > 2)
    error('input data has too many columns!');
end

% Assign information
time_series.data = in_data;
time_series.name = ts_name;
time_series.short_name = ts_short_name;
time_series.unit = ts_unit;

% Save structure to file
if ismac
    save_dir = 'Timeseries_datasets/';
elseif ispc
    save_dir = 'Timeseries_datasets\';
end
out_name = strcat(save_dir,data_name,'.mat');
    
save(out_name,'time_series');

disp(['Created ' out_name])

