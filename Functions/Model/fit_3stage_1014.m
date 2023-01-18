%
% out = fit_3stage_1014(sample_data,model_time_bnds,bur_dur_bnds,n_iterations)
% out = fit_3stage_1014(sample_data,model_time_bnds,bur_dur_bnds,n_iterations,model_interval,plot_fig,burial_frac_bnds)
%
% Computes the predicted nuclide concentration for each sample and nuclide,
% based on the given time-series driver and floating parameters. Then 
% calculates the misfit with measured nuclide concentrations.
%
% sample_data is a required struct, created using get_data_1014.m.
%
% model_time_bnds should be the bounds of the total model time (years 
% before present) [min max], used for 3-stage driver.
%
% bur_dur_bnds should be the bounds of the burial duration [min max], used 
% for 3-stage driver.
%
% n_iterations is the number of model iterations.
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
% [min max], used for the 3-stage driver.
%
% Output are the misfits between predicted and measured nuclide
% concentrations, for each sample.
%
%
%%

function out = fit_3stage_1014(sample_data,model_time_bnds,bur_dur_bnds,n_iterations,model_interval,plot_fig,burial_frac_bnds)

  % Check inputs
  if (nargin < 4 || nargin > 7)
      error('fit_3stage_1014 has wrong number of inputs!');
  end
  if (nargin < 5) || isempty(model_interval)
      model_interval = 1000;
  end
  if (model_interval ~= 1000 && model_interval ~= 100 && model_interval ~= 10)
      error('model_interval must be 10, 100 or 1000 years!');
  end
  if (nargin < 6)
      plot_fig = [];
  end
  if (nargin < 7)
      burial_frac_bnds = [];
  end
  
  if model_interval == 1000
      r_fac = -3;
  elseif model_interval == 100
      r_fac = -2;
  elseif model_interval == 10
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
          data{a}.meanage14 = sample_data.ages.C14(a,1);
      end
      ages14(a)=data{a}.meanage14;
  end
  
  
  % Determine the random model paraemters
  rand_model_time = param_rand_gen(model_time_bnds(1),model_time_bnds(2),n_iterations);
  rand_bur_dur = param_rand_gen(bur_dur_bnds(1),bur_dur_bnds(2),n_iterations);
  
  
  % Calculate for each scenario
  disp('Finding best fit scenario...');
  
  for b = 1:n_iterations
      
      bur_dur = rand_bur_dur(b);
      model_time = rand_model_time(b);
      
      
      % Generate exposure and burial periods
      recent_exp = mean(ages14);
      exposed_or_not = three_stage_logical(recent_exp,bur_dur,model_time,model_interval);
      
          
      if ~isempty(burial_frac_bnds)
%           % Determine burial fractions for each sample depending on
%           % elevation order
%           burial_frac_int = 0.05; % Interval for burial fractions within given range
%           burial_fracs = burial_frac_bnds(1):burial_frac_int:burial_frac_bnds(2)+burial_frac_int; % Burial ractions
%           
%           [~,pos_idx] = sort(sample_data.position,'descend'); % Get sample position order
%           sample_arr = 1:length(sample_data.position); sample_arr_rev = fliplr(sample_arr);
%           for c = sample_arr % Create array for each sample, staggered by position
%               this_pos_idx = pos_idx(c);
%               sample_burial_fracs(:,this_pos_idx) = burial_fracs(sample_arr(c):end-sample_arr_rev(c));
%           end
          
          % Determine burial fractions for each sample depending on
          % relative elevation difference
          n_bur_frac = 10; % Set number of iterations to test burial fraction
          bur_frac_arr = linspace(burial_frac_bnds(1),burial_frac_bnds(2),n_bur_frac);
          pos_scaled = rescale(sample_data.position);
          for c = 1:length(bur_frac_arr)
              sample_burial_fracs(c,:) = bur_frac_arr(c)*pos_scaled;
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
              
              for d = 1:length(sample_burial_fracs(:,1))
                  
                  this_frac_buried = sample_burial_fracs(d,c);
                  this_recent_exp = data_s.meanage14;
                  
                  % Generate exposure and burial periods
                  this_exposed_or_not = three_stage_logical(this_recent_exp,bur_dur,model_time,model_interval,this_frac_buried);
                  
                  % Calculate predicted concentrations and misfit
                  [this_burfrac_predN,this_s_misfit] = fit_sample(this_exposed_or_not,data_s);
                  
                  burfrac_predN{d} = this_burfrac_predN;
                  this_burfrac_misfit(d) = this_s_misfit;
                  
              end
               
              % Get the bestfit burial fraction
              if numel(unique(this_burfrac_misfit))==numel(this_burfrac_misfit)
                  burfrac_min_misfit = min(this_burfrac_misfit);
                  this_s_predN = burfrac_predN(burfrac_min_misfit==this_burfrac_misfit);
                  this_s_bestfit_burfrac = sample_burial_fracs((burfrac_min_misfit==this_burfrac_misfit),c);
              else % Use the first result if no unique bestfit
                  burfrac_min_misfit = min(this_burfrac_misfit);
                  this_s_predN = burfrac_predN(1);
                  this_s_bestfit_burfrac = sample_burial_fracs(1,c);
              end
              this_scenario_predN{c} = this_s_predN{1};
              this_scenario_misfit(c) = burfrac_min_misfit;
              this_scenario_burfrac(c) = this_s_bestfit_burfrac;
          end
      end
      if isempty(burial_frac_bnds)
          scenario_misfits(b) = mean(this_scenario_misfit,'omitnan');
      else
          if numel(unique(this_burfrac_misfit))==numel(this_burfrac_misfit)
              scenario_misfits(b) = mean(this_scenario_misfit,'omitnan');
              scenario_bestfit_burfracs(b,:) = this_scenario_burfrac;
          else
              scenario_misfits(b) = NaN;
              scenario_bestfit_burfracs(b,:) = NaN;
              continue
          end
      end
      
      % Plot concentrations for scenario
      if ~isempty(plot_fig)
          if b>1
              delete(plot_fig(1,2).Children(1:length(this_scenario_predN)));
          end
          plot_Scenario(exposed_or_not,sample_data,this_scenario_predN,plot_fig)
      end
      
  end
  
  
  % Determine bestfit scenario
  bestfit_log = find(min(scenario_misfits) == scenario_misfits);
  bestfit_bur_dur = round(rand_bur_dur(bestfit_log),r_fac);
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
  disp(['bur_dur (years)  ' int2str(bestfit_bur_dur(1)) ]);
  disp(['model_time (years)  ' int2str(bestfit_model_time(1)) ]);
  
  % Export
  out.misfit_min = min(scenario_misfits);
  out.bestfit_bur_dur = bestfit_bur_dur;
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
        ts_plot = plot_expbur_time(exposed_or_not,plot_fig);
        
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
