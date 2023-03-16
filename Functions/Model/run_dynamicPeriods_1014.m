%
% out = run_dynamicPeriods_1014(sample_data,expo_mids,modeltime_initial,expodur_initial,misfit_type)
% out = run_dynamicPeriods_1014(sample_data,expo_mids,modeltime_initial,expodur_initial,misfit_type,model_interval,plot_fig,burialfrac_bnds)
%
% Computes the predicted nuclide concentration for each sample and nuclide,
% based on a multi-stage exposure-burial history, using floating
% parameters with an optimised solver. Misfit is calculated using the
% measured nuclide concentrations.
%
% sample_data is a required struct, created using get_data_1014.m.
%
% expo_mids is a required vector of mid-points for the exposure periods (in
% kilo years ago).
%
% modeltime_initial should be a starting value of the total model time 
% (years before present), used for dynamic periods driver.
%
% expodur_initial should be the starting value of the burial duration, used 
% for dynamic periods driver.
%
% misfit_type specifies the method of deriving misfit ('all','minmax',
% 'minBe','maxBe','minC','maxC'). For example, 'all' uses mean of all 
% samples, 'minmax' uses the mean of the min and max sample misfits, and 
% 'maxBe' to use the misfit from the sample with the maximum 10-Be 
% concentration.
%
% model_interval is an optional input, used in calculations - 10, 100 or 
% 1000 years (default is 1000).
%
% plot_fig is an optional input of the figure handle generated with
% plot_concs_1014.m. If included, the time series, exposure/burial periods
% and corresponding predicted nuclide concentrations are plotted for each
% modelled scenario.
%
% burial_frac_bnds is an optional input, specifying the fraction of 
% exposure period that a sample could be buried (multiples of 0.05) 
% [min max], used for the dynamic periods driver.
%
% Output is the bestfit model time and exposure duration, and the 
% corresponding misfit between predicted and measured nuclide concentrations.
%
%
%%

function out = run_dynamicPeriods_1014(sample_data,expo_mids,modeltime_initial,expodur_initial,misfit_type,model_interval,plot_fig,burialfrac_bnds)

  % Check inputs
  if (nargin < 5 || nargin > 8)
      error('run_dynamicPeriods_1014 has wrong number of inputs!');
  end
  if (nargin < 6) || isempty(model_interval)
      model_interval = 1000;
  end
  if (model_interval ~= 1000 && model_interval ~= 100 && model_interval ~= 10)
      error('model_interval must be 10, 100 or 1000 years!');
  end
  if (nargin < 7)
      plot_fig = [];
  end
  if (nargin < 8)
      burialfrac_bnds = [];
  end

  
  % Re-organise data for each sample
  sample_data.logical_1014 = any(sample_data.logical_10 & sample_data.logical_14); 
  for a = 1:length(sample_data.s)
      data{a}.name = sample_data.s{a}.name{1};
      data{a}.cover = sample_data.cover;
      %data{a}.ee = ee;
      data{a}.pp = sample_data.pp;
      data{a}.sf10 = sample_data.sf1026(a);
      data{a}.cp10 = sample_data.cp1026(a);
      data{a}.sf14 = sample_data.sf14(a);
      data{a}.cp14 = sample_data.cp14(a);
      data{a}.top_z_gcm2 = sample_data.s{a}.top_z_gcm2;
      data{a}.bottom_z_gcm2 = sample_data.s{a}.bottom_z_gcm2;
      data{a}.weight10 = sample_data.s{a}.weight10;
      data{a}.weight14 = sample_data.s{a}.weight14;
      data{a}.nuclide10 = sample_data.s{a}.nuclide10;
      data{a}.nuclide14 = sample_data.s{a}.nuclide14;
      if sample_data.s{a}.nuclide10 == 1
          data{a}.N10 = sample_data.s{a}.N10;
          data{a}.dN10 = sample_data.s{a}.dN10;
      end
      if sample_data.s{a}.nuclide14 == 1
          data{a}.N14 = sample_data.s{a}.N14;
          data{a}.dN14 = sample_data.s{a}.dN14;
          data{a}.meanage14 = sample_data.ages.C14(a,1);
      end
      sample_data.ages14(a)=data{a}.meanage14;
  end
  sample_data.org_data = data;
  
  
  % Set optimisation rules
  opts = optimset('fminsearch');
  opts = optimset(opts,'TolFun',1,'TolX',1);%,'Display','iter');
  
  
  % Find bestfit scenario and export result
  disp('Finding best fit scenario...');
  
  [optX,fmin] = fminsearch(@(X) fit_opt_dynamicPeriods_1014(X,sample_data,expo_mids,model_interval,misfit_type,plot_fig,burialfrac_bnds),[modeltime_initial,expodur_initial],opts);
  
  disp('');
  disp('Best fit scenario:');
  disp(['expo_dur (years)  ' int2str(optX(2)) ]);
  disp(['model_time (years)  ' int2str(optX(1)) ]);
  
  
  % Export
  out.misfit_min = fmin;
  out.bestfit_expo_dur = optX(2);
  out.bestfit_model_time = optX(1);
  
  % Calculate reduced chi-squared from bestfit
  if length(sample_data.s) > 2
      DOF = length(sample_data.s)-2;
      reduced_chi2 = fmin ./ DOF;

      disp(['Reduced chi-squared is ' sprintf('%0.2f',reduced_chi2) ' for ' int2str(DOF) ' DOF']);
      out.bestfit_reduced_chi2 = reduced_chi2;
  end
  
end
