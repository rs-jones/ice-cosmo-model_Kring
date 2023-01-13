%
% out = fit_timeseries_1014(sample_data,driver,threshold_frac_bnds,model_time_bnds,n_iterations)
% out = fit_timeseries_1014(sample_data,driver,threshold_frac_bnds,model_time_bnds,n_iterations,plot_fig)
%
% Computes the predicted nuclide concentration for each sample and nuclide,
% based on the given time-series driver and floating parameters. Then 
% calculates the misfit with measured nuclide concentrations.
%
% sample_data is a required struct, created using get_data_1014.m.
%
% driver is a required struct containing data necessary for computing 
% periods of exposure/burial (time_series_logical).
%
% threshold_frac_bnds should be the bounds of the threshold fraction 
% [min max], used for time series driver.
%
% model_time_bnds should be the bounds of the total model time (years 
% before present) [min max], used for time series driver.
%
% n_iterations is the number of model iterations.
%
% plot_fig is an optional input of the figure handle generated with
% plot_concs_1014.m. If included, the time series, exposure/burial periods
% and corresponding predicted nuclide concentrations are plotted for each
% modelled scenario.
%
% burial_frac_bnds is an optional input, specifying the fraction of 
% exposure period that a sample could be buried (multiples of 0.05) 
% [min max], used for the time series driver.
%
% Output are the misfits between predicted and measured nuclide
% concentrations, for each sample.
%
%
%%

function out = fit_timeseries_1014(sample_data,driver,threshold_frac_bnds,model_time_bnds,n_iterations,plot_fig,burial_frac_bnds)

  % Check inputs
  if (nargin < 5 || nargin > 7)
      error('fit_timeseries_1014 has wrong number of inputs!');
  end
  if (nargin < 6)
      plot_fig = [];
  end
  if (nargin < 7)
      burial_frac_bnds = [];
  end
  
  if driver.model_interval == 1000
      r_fac = -3;
  elseif driver.model_interval == 100
      r_fac = -2;
  elseif driver.model_interval == 10
      r_fac = -1;
  end
  

  % Determine whether samples have multiple nuclide measurements
  logical_1014 = any(sample_data.logical_10 & sample_data.logical_14); 
  
  
  % Re-organise data for each sample
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
      end
  end
  
  
  % Determine the random model paraemters
  rand_threshold_frac = param_rand_gen(threshold_frac_bnds(1),threshold_frac_bnds(2),n_iterations);
  rand_model_time = param_rand_gen(model_time_bnds(1),model_time_bnds(2),n_iterations);

  
  % Calculate for each scenario
  disp('Finding best fit scenario...');
  
  for b = 1:n_iterations
      
      threshold_frac = rand_threshold_frac(b);
      model_time = rand_model_time(b);
      time_series = driver.time_series;
      data_use = driver.data_use;
      model_interval = driver.model_interval;
      
      
      % Generate exposure and burial periods
      exposed_or_not = time_series_logical(time_series,threshold_frac,model_time,data_use,model_interval);
      
      
      % Determine burial fractions for each sample depending on
      % relative position
      if ~isempty(burial_frac_bnds)
          burial_frac_int = 0.05; % Interval for burial fractions within given range
          burial_fracs = burial_frac_bnds(1):burial_frac_int:burial_frac_bnds(2)+burial_frac_int; % Burial ractions
          
          [~,pos_idx] = sort(sample_data.position,'descend'); % Get sample position order
          sample_arr = 1:length(sample_data.position); sample_arr_rev = fliplr(sample_arr);
          for c = sample_arr % Create array for each sample, staggered by position
              this_pos_idx = pos_idx(c);
              sample_burial_fracs(:,this_pos_idx) = burial_fracs(sample_arr(c):end-sample_arr_rev(c));
          end
          
      end
      
      
      % Evaluate each sample
      for c = 1:length(data)
          
          data_s = data{c};
          
          
          % If no burial fractions are specified, run once per sample
          if isempty(burial_frac_bnds)
              
              % Calculate predicted concentrations and misfit
              [this_s_predN,this_s_misfit] = fit_sample(exposed_or_not,data_s);
              
              this_scenario_predN{c} = this_s_predN;
              this_scenario_misfit(c) = this_s_misfit;
          
          % Otherwise, do range of burial fractions for each sample
          else
              
              for d = 1:length(sample_burial_fracs)
                  
                  this_frac_buried = sample_burial_fracs(d,c);
                  
                  % Generate exposure and burial periods
                  this_exposed_or_not = time_series_logical(time_series,threshold_frac,model_time,data_use,model_interval,this_frac_buried);
                  
                  % Calculate predicted concentrations and misfit
                  [this_burfrac_predN,this_s_misfit] = fit_sample(this_exposed_or_not,data_s);
                  
                  burfrac_predN{d} = this_burfrac_predN;
                  this_burfrac_misfit(d) = this_s_misfit;
                  
              end
                            
              if numel(unique(this_burfrac_misfit))==numel(this_burfrac_misfit)
                  burfrac_min_misfit = min(this_burfrac_misfit);
                  this_s_predN = burfrac_predN(burfrac_min_misfit==this_burfrac_misfit);
                  this_s_bestfit_burfrac = sample_burial_fracs((burfrac_min_misfit==this_burfrac_misfit),c);
                  
                  this_scenario_predN{c} = this_s_predN{1};
                  this_scenario_misfit(c) = burfrac_min_misfit;
                  this_scenario_burfrac(c) = this_s_bestfit_burfrac;
              end
          end
          
      end
      scenario_misfits(b) = mean(this_scenario_misfit,'omitnan');
      if ~isempty(burial_frac_bnds)
          scenario_bestfit_burfracs(b,:) = this_scenario_burfrac;
      end
      
      % Plot concentrations for scenario
      if ~isempty(plot_fig)
          if b>1
              delete(plot_fig(1,2).Children(1:length(this_scenario_predN)));
          end
          plot_Scenario(exposed_or_not,sample_data,this_scenario_predN,plot_fig)
      end
      
%       disp(['Mean misfit of  ' int2str(mean(this_scenario_misfit)) '  from scenario:']);
%       disp(['threshold_frac  ' sprintf('%0.2f',threshold_frac) ]);
%       disp(['model_time (years)  ' int2str(model_time) ]);
      
  end
  
  
  % Determine bestfit scenario
  bestfit_log = find(min(scenario_misfits) == scenario_misfits);
  bestfit_threshold_frac = rand_threshold_frac(bestfit_log);
  bestfit_model_time = round(rand_model_time(bestfit_log),r_fac);
  if ~isempty(burial_frac_bnds)
      bestfit_burfracs = scenario_bestfit_burfracs(bestfit_log,:);
  end
  
  % Calculate reduced chi-squared from bestfit
  if length(sample_data.s) > 2
      DOF = length(sample_data.s)-2;
      reduced_chi2 = min(scenario_misfits) ./ DOF;
  end
  
  disp('');
  disp('Best fit scenario:');
  disp(['threshold_frac  ' sprintf('%0.2f',bestfit_threshold_frac(1)) ]);
  disp(['model_time (years)  ' int2str(bestfit_model_time(1)) ]);
  
  % Export
  out.misfit_min = min(scenario_misfits);
  out.bestfit_threshold_frac = bestfit_threshold_frac;
  out.bestfit_model_time = bestfit_model_time;
  if ~isempty(burial_frac_bnds)
      out.bestfit_burfracs = bestfit_burfracs;
  end
  
  if length(sample_data.s) > 2
      disp(['Reduced chi-squared is ' sprintf('%0.2f',reduced_chi2) ' for ' int2str(DOF) ' DOF']);
      out.bestfit_reduced_chi2 = reduced_chi2;
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%% Evaluate each sample %%%%%%%%%%%%%%%%%%%%%%%%%
  
    function [predN,misfit] = fit_sample(exposed_or_not,data_s)
        
        ee = [];
        cover = [];
        
        % Calculate predicted concentrations
        predN = forward_model_1014(exposed_or_not,data_s,ee,cover);
        
        
        % Calculate misfit
        if logical_1014
            
            misfit_10 = ((predN.sum.N10 - data_s.N10)./data_s.dN10) .^2;
            misfit_14 = ((predN.sum.N14 - data_s.N14)./data_s.dN14) .^2;
            misfit = mean([misfit_10,misfit_14]); % Average misfits of each nuclide
        else
            if data_s.nuclide10 == 1
                misfit = ((predN.sum.N10 - data_s.N10)./data_s.dN10) .^2;
            end
            if data_s.nuclide14 == 1
                misfit = ((predN.sum.N14 - data_s.N14)./data_s.dN14) .^2;
            end
        end
        
    end
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot scenario %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    function plot_Scenario(exposed_or_not,data,scenario_predN,plot_fig)
        
        figure(plot_fig(1));
        
        % Plot time-series data
        ts_plot = plot_time_series(exposed_or_not,plot_fig);
        
        % Plot predicted concentrations for samples
        for d = 1:length(scenario_predN)
            this_predN = scenario_predN{d};
            plot_concs_1014_pred(data,plot_fig,this_predN,[]);
        end
        
        % Plot predicted concentration pathways for samples
        %plot_pred_pathway(sample_data,plot_fig,predN,[]); % FIX
        
        % Plot exposure/burial periods
        plot_ExpBur(exposed_or_not,ts_plot,plot_fig);
        
        drawnow;
        
    end


end
