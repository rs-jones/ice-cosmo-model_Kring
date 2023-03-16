%% ANALYSE MULTIPLE NUCLIDE CONCENTRATIONS FROM SURFACE SAMPLES USING A TIME-SERIES DRIVER
%
% Finds the best fit exposure/burial scenario to explain the measured
% nuclide concentations. Plots the corresponding history, and predicted
% concentrations on a two-isotope diagram.
%
% Uses time-series data, with floating model time (as years before present)
% and threshold (as a fraction of the maximum time series value), to 
% compute the exposure/burial history.
%
% To be used with nuclides measured in surface samples (e.g. bedrock, 
% boulders).
%
% Input sample data must include the following information:
% 1. Sample name
% 2. Latitude (decimal degrees)
% 3. Longitude (decimal degrees)
% 4. Elevation (m asl)
% 5. Pressure (hPa) (zero if not known)
% 6. Relative position (distance from terminus, km; elevation above ice, m)
% 7. Sample thickness (cm)
% 8. Bulk density (g/cm^3)
% 9. Shielding factor for terrain, snow, etc. (unitless)
% 10. Sample 10-Be concentration (atoms of 10-Be/g)
% 11. Sample 10-Be concentration 1 sigma uncertainty (atoms of 10-Be/g)
% 12. Sample 14-C concentration (atoms of 14-C/g)
% 13. Sample 14-C concentration 1 sigma uncertainty (atoms of 14-C/g)
% 14. Top depth of sample (cm)
% 15. Bottom depth of sample (cm)
% 16. 10-Be final mineral weight (g)
% 17. 14-C final mineral weight (g)
% 18. Year the sample was collected (calendar year)
%
% Optional data should include the following information:
% 19. Sample 10-Be exposure age (mean; years)
% 20. Sample 10-Be exposure 1 sigma uncertainty (internal; years)
% 21. Sample 10-Be exposure 1 sigma uncertainty (external; years)
% 22. Sample 14-C exposure age (mean; years)
% 23. Sample 14-C exposure 1 sigma uncertainty (internal; years)
% 24. Sample 14-C exposure 1 sigma uncertainty (external; years)
% 25. Scaling model used (i.e. 'DE','DU','LI','ST','LM','LSD'/'SA','LSDn'/'SF')
%
%
%%
clear % Start fresh
addpath(genpath(pwd))

% SET inputs
input_name = 'MtKring_BurExp_input'; % Name used for for sample data .xlsx
scaling_model = 'ST'; % Scaling model - time-independent scaling (e.g. Stone) only currently available


% Load and assign time-series data
load('SST_SO_stack.mat'); % Specify forcing time series (create using Create_timeseries.m)
%load('EDC_atm_temp.mat'); % Example time series (Epica Dome C atm. temp.)
%load('LR04_stack.mat'); % Example time series (Lisiecki & Raymo '04 benthic stack)
%load('PD09_ice_vol.mat'); % Example time series (Pollard & DeConto '09 Antarctic ice volume)

data_use = 'normal'; % OPTIONALLY set data use (default is 'normal' if []). 
    % If 'opposite' is set, exposure intervals become burial intervals, and vice versa. 
model_interval = 1000; % OPTIONALLY set model interval (10, 100 or 1000 years; 
    % default is 1000 years), based on the resolution of the time-series data. 


%% Get Data

% Load sample data
sample_data = get_data_1014(input_name);

% Get parameters for samples
sample_data = get_pars_1014(sample_data,scaling_model);

% Assign time series data to model driver
driver.time_series = time_series;
driver.data_use = data_use;
driver.model_interval = model_interval;


%% Determine Best Fit Scenarios

% SET model bounds
threshold_frac_bnds = [0.2 0.8];    % Threshold fraction for time series driver (lower and upper)
model_time_bnds = [20000 100000];   % Total model time (lower and upper)
burial_frac_bnds = [0 .9];          % Fraction of exposure period that a sample could be buried (lower and upper)
n_iterations = 1000;                % Number of model iterations


% Generate figure of sample concentrations and plot handles
expo_intervals = [25,75,50,100]; % Specify exposure intervals to show on plot
bur_intervals = [10,20]; % Specify burial intervals to show on plot
x_lim = [0,1.2e5]; % Specify 10Be axis limits
y_lim = [0.012,1]; % Specify 14C/10Be axis limits
add_names = 0; % Add sample names to plot? (0=no, 1=yes)

close all; % Close existing figure, if open
fig_h = plot_concs_1014(sample_data,2,1,expo_intervals,bur_intervals,x_lim,y_lim,add_names);


% Find best fit for each sample
best_fits = run_timeseries_1014(sample_data,driver,threshold_frac_bnds,model_time_bnds,n_iterations,fig_h,burial_frac_bnds);


% Plot bestfit concentrations
close;
fig_bf = plot_bestfit_ts_1014(sample_data,driver,best_fits,expo_intervals,bur_intervals,x_lim,y_lim);


% Export bestfit results
save_name = 'MtKring_SST_n1000'; % Specify file name
%save(strcat(save_name,'_output.mat'),'best_fits'); % Save data
saveas(fig_bf,strcat(save_name,'.png')); % Save figure

