%% ANALYSE MULTIPLE NUCLIDE CONCENTRATIONS FROM SURFACE SAMPLES USING A THREE-STAGE HISTORY
%
% Finds the best fit exposure/burial scenario to explain the measured
% nuclide concentations. Plots the corresponding history, and predicted
% concentrations on a two-isotope diagram.
%
% Uses the floating parameters of total exposure-burial time and time 
% buried (in years), to compute a three-stage (exposure-burial-exposure) 
% history.
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
% 19. Sample 10-Be exposure age - Optional (can use NaN)
% 20. Sample 10-Be exposure 1 sigma uncertainty - Optional (can use NaN)
% 21. Sample 10-Be exposure 1 sigma uncertainty - Optional (can use NaN)
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


% Load sample data
sample_data = get_data_1014(input_name);

% Get parameters for samples
sample_data = get_pars_1014(sample_data,scaling_model);


%% Determine Best Fit Scenarios

% SET model inputs
modeltime_initial = 100000;        % Total model time (years, starting value)
burdur_initial = 50000;            % Duration of burial (years; starting value)
burialfrac_bnds = [0 .9];          % Fraction of exposure period that a sample could be buried (lower and upper)
model_interval = [];               % Optionally set model interval used in calculations - 10, 100 or 1000 years (default is 1000)
misfit_type = 'all';               % Set method of deriving misfit - 'all','minmax','minBe','maxBe','minC','maxC'
                                       % 'all' to use mean of all samples
                                       % 'minmax' to use the mean of the min and max sample misfits
                                       % 'maxBe' to use the misfit from the sample with the maximum 10-Be concentration

                                         
% Generate figure of sample concentrations and plot handles
expo_intervals = [25,75,50,100]; % Specify exposure intervals to show on plot
bur_intervals = [10,20]; % Specify burial intervals to show on plot
x_lim = [0,1.2e5]; % Specify 10Be axis limits
y_lim = [0.012,1]; % Specify 14C/10Be axis limits
add_names = 0; % Add sample names to plot? (0=no, 1=yes)

close all; % Close existing figure, if open
fig_h = plot_concs_1014(sample_data,2,1,expo_intervals,bur_intervals,x_lim,y_lim,add_names);


% Find best fit for each sample
best_fits = run_3stage_1014(sample_data,modeltime_initial,burdur_initial,misfit_type,model_interval,fig_h,burialfrac_bnds);


% Export bestfit results
save_name = 'MtKring_3stage_bestfit'; % Specify file name
%save(strcat(save_name,'_output.mat'),'best_fits'); % Save data
saveas(gcf,strcat(save_name,'.png'),'png'); % Save figure

