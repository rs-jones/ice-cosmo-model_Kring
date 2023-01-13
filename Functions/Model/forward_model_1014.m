%
% predN = forward_model_1014(exposed_or_not,data_s)
% predN = forward_model_1014(exposed_or_not,data_s,ee,cover)
%
% Calculates predicted nuclide concentrations for each nuclide
% by integrating over time and depth, for a single sample.
%
% It uses generated time intervals of exposure/burial (exposed_or_not).
%
% data_s is a required struct of data for a single sample, created in 
% fit_timeseries_1014.m.
%
% ee is an optional struct of erosion rates, for periods of exposure
% (ee.e_expo) and burial (ee.e_bur). Defaults are zero.
%
% cover is an optional struct containing surface cover depth (cover.z).
% Default is zero.
%
% Output is a struct of predicted nuclide concentrations for all time 
% intervals and summed intervals, and then for each nuclide.
%
% THIS CURRENTLY ONLY WORKS FOR TIME-INDEPENDENT SCALING SCHEMES, AND
% EXCLUDES EROSION DURING BURIAL PERIODS
%
%
%%

function predN = forward_model_1014(exposed_or_not,data_s,ee,cover)

  %warning('off','MATLAB:integral2:funVectorization');

  % Check inputs
  if (nargin < 2 || nargin > 4)
      error('forward_model_1014 has wrong number of inputs!');
  end
  
  if (nargin < 3) || isempty(ee)
      ee.e_expo = 0; % If ee is missing, set erosion rate during exposure to zero
      ee.e_bur = 0;  % If ee is missing, set erosion rate during burial to zero
  end
  if (nargin < 3) || isempty(cover)
      try 
          cover.z = data_s.cover.z; % If cover is missing, find in sample_data
      catch
          cover.z = 0; % If cover cannot be found in sample_data, set cover depth to zero
      end
  end

  
  if 1; % Option (1 = on) -- CHANGE SO DEPENDENT ON WHETHER TIME-DEPENDENT SCALING IS SET
      interval_time = exposed_or_not.reduced_interval_time;
      logical = exposed_or_not.reduced_logical;
      cum_t_expo = exposed_or_not.cum_t_expo;
      timeBP = exposed_or_not.reduced_time;
  else
      interval_time = exposed_or_not.interval_time;
      logical = exposed_or_not.logical;
      cum_t_expo = interval_time;
      timeBP = exposed_or_not.time;
  end
  

  n_intervals = numel(logical);
  
  N_10 = zeros(n_intervals,1); % row for each time interval
  N_14 = N_10;
  
  % Calculate for each interval
  for a = 1:n_intervals
      
      % Exposure interval
      if logical(a) == 1
          
          
          % Compute depths at the end of this exposure period
          if a < n_intervals
              this_top_z = data_s.top_z_gcm2 + ee.e_expo.*cum_t_expo(a+1) + cover.z;
              this_bottom_z = data_s.bottom_z_gcm2 + ee.e_expo.*cum_t_expo(a+1) + cover.z;
          else
              this_top_z = data_s.top_z_gcm2;
              this_bottom_z = data_s.bottom_z_gcm2;
          end
         
          if data_s.nuclide10 == 1
              
              % Define production rate functions
              l = data_s.pp.lambda10Be;
              pr_func = @(z,t) PR_Z((z + ee.e_expo.*t),data_s.pp,data_s.sf10{1},data_s.cp10{1},10) .* exp(-l.*t);
              
              % Integrate over time and depth for each sample
              Integrated_N = zeros(size(this_top_z));
              for b = 1:length(this_top_z)
                  % Integrate each sample from top depth to bottom depth and from zero to the exposure time in this time step.
                  % Divide integral by sample thickness.
                  Integrated_N(b) = integral2(pr_func,this_top_z(b),this_bottom_z(b),0,interval_time(a),'RelTol',1e-3,'AbsTol',1e-3) ./ (this_bottom_z(b) - this_top_z(b));
              end           
              
              % Average by mineral weight and account for subsequent burial
              this_weight = data_s.weight10;
              if a < n_intervals
                  N_10(a) =  ( (sum(Integrated_N.*this_weight))./sum(this_weight) ) .* exp(-l.*timeBP(a+1));
              else
                  N_10(a) =  ( (sum(Integrated_N.*this_weight))./sum(this_weight) );
              end
          end
          
          if data_s.nuclide14 == 1
              
              % Define production rate functions
              l = data_s.pp.lambda14C;
              pr_func = @(z,t) PR_Z((z + ee.e_expo.*t),data_s.pp,data_s.sf14{1},data_s.cp14{1},14) .* exp(-l.*t);
              
              % Integrate over time and depth for each sample
              Integrated_N = zeros(size(this_top_z));
              for b = 1:length(this_top_z)
                  % Integrate each sample from top depth to bottom depth and from zero to the exposure time in this time step.
                  % Divide integral by sample thickness.
                  Integrated_N(b) = integral2(pr_func,this_top_z(b),this_bottom_z(b),0,interval_time(a),'RelTol',1e-3,'AbsTol',1e-3) ./ (this_bottom_z(b) - this_top_z(b));
              end
              
              % Average by mineral weight and account for subsequent burial
              this_weight = data_s.weight14;
              if a < n_intervals
                  N_14(a) =  ( (sum(Integrated_N.*this_weight))./sum(this_weight) ) .* exp(-l.*timeBP(a+1));
              else
                  N_14(a) =  ( (sum(Integrated_N.*this_weight))./sum(this_weight) );
              end
          end
          
      else
          % Burial period - do nothing
          
      end
      
  end

  
  predN.all.N10 = N_10;
  predN.all.N14 = N_14;
    
  % Sum all intervals
  predN.sum.N10 = sum(N_10);
  predN.sum.N14 = sum(N_14);


end
