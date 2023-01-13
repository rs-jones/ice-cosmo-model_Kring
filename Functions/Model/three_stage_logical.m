%
% exposed_or_not = three_stage_logical(recent_exp,bur_dur,model_time)
% exposed_or_not = three_stage_logical(recent_exp,bur_dur,model_time,model_interval,frac_buried)
%
% Creates a logical of time intervals that are either exposed or buried
% based on a simple three-stage exposure-burial-exposre history.
%
% Output is to be used to model nuclide concentrations through time.
%
% recent_exp (in years) is the duration of recent exposure.
%
% bur_dur (in years) is the duration of preceding continuous burial.
%
% model_time is the total time of exposure-burial-exposure (in years).
%
% model_interval (in years) is optional. Default is 1000 years.
%
% frac_buried is the fraction of an exposure period that a sample is buried
% (optional).
%
% The computed logical, and reduced times and logical are exported to a 
% single structure.
%
%
%%

function exposed_or_not = three_stage_logical(recent_exp,bur_dur,model_time,model_interval,frac_buried)

  
  % Check inputs
  if (nargin < 3 || nargin > 5)
      error('three_stage_logical has the wrong number of inputs!');
  end
  
  if (nargin < 3) || isempty(model_interval)
      model_interval = 1000;
  end
  if (model_interval ~= 1000 && model_interval ~= 100 && model_interval ~= 10)
      error('model_interval must be 10, 100 or 1000 years!');
  end
  if (nargin < 5) || isempty(frac_buried)
      frac_buried = 0;
  end
  
  if model_interval == 1000
      r_fac = -3;
  elseif model_interval == 100
      r_fac = -2;
  elseif model_interval == 10
      r_fac = -1;
  end
  
  
  EBE_time = round(model_time,r_fac);
  recent_exp_dur = round(recent_exp,r_fac);
  buried_dur = round(bur_dur,r_fac);
  
  % Determine initial exposure time
  init_exp_dur = EBE_time-buried_dur-recent_exp_dur;
  
  if init_exp_dur<1
      init_exp_dur = 1;
      EBE_time = init_exp_dur+buried_dur+recent_exp_dur;
  end
  time = fliplr(0:model_interval:EBE_time);
  
  
  % Compute logical
  init_exposure = ones(1,init_exp_dur/model_interval);
  burial = zeros(1,buried_dur/model_interval);
  recent_exposure = ones(1,recent_exp_dur/model_interval);
  exposed_or_not.logical = logical([init_exposure burial recent_exposure 1]);
  
  
  % Append time and interval
  exposed_or_not.time = time;
  exposed_or_not.values = [init_exposure burial recent_exposure 1];
  exposed_or_not.interval_time = zeros(1,numel(exposed_or_not.logical)) + model_interval;
  
  
  % Account for fraction of exposure period that is buried, and
  % reduce exposed_or_not to single intervals of alternate exposure/burial
  if frac_buried>0
      expbur_change = find(diff([-1 exposed_or_not.logical -1]) ~= 0); % Find points along logical at either end of exposure periods
      interval_time = diff(expbur_change) .* model_interval; % Calculate time difference between these points
      expbur_change = expbur_change(1:end-1);
      n_intervals = numel(interval_time);
      for ii = 1:n_intervals
          reduced_time(ii) = exposed_or_not.time(expbur_change(ii)); % Determine new times
          reduced_logical(ii) = (exposed_or_not.logical(expbur_change(ii)) == 1); % Determine new logical
      end
      exposed_or_not.reduced_time_zeroBurFrac = reduced_time;
      exposed_or_not.reduced_logical_zeroBurFrac = reduced_logical;
      exp_times_buried = round(interval_time * frac_buried); % Calculate times buried for each interval
      exp_times_buried(~reduced_logical) = 0; % Keep only exposure interval times
      for iii = 1:n_intervals
          if reduced_logical(iii)
              interval_time(iii) = interval_time(iii) - exp_times_buried(iii); % Remove burial time from exposure period
              if iii<n_intervals
                  interval_time(iii+1) = interval_time(iii+1) + exp_times_buried(iii);  % Add burial time to burial period
                  reduced_time(iii+1) = reduced_time(iii+1) + exp_times_buried(iii); % Adjust start of burial period
              end
          end
          reduced_values(iii) = exposed_or_not.values(expbur_change(iii)); % Determine new values
      end
      t_expo = interval_time; t_expo(~reduced_logical) = 0; % Convert to time exposed (make burial periods zero)
      cum_t_expo = fliplr(cumsum(fliplr(t_expo))); % Calculate cumulative exposure time (from present)
  else
      % Reduce exposed_or_not to single intervals of alternate exposure/burial for improved efficiency
      expbur_change = find(diff([-1 exposed_or_not.logical -1]) ~= 0); % Find points along logical at either end of exposure periods
      interval_time = diff(expbur_change) .* model_interval; % Calculate time difference between these points
      expbur_change = expbur_change(1:end-1);
      n_intervals = numel(interval_time);
      for ii = 1:n_intervals
          reduced_time(ii) = exposed_or_not.time(expbur_change(ii)); % Determine new times
          reduced_logical(ii) = (exposed_or_not.logical(expbur_change(ii)) == 1); % Determine new logical
          reduced_values(ii) = exposed_or_not.values(expbur_change(ii)); % Determine new values
      end
      t_expo = interval_time; t_expo(~reduced_logical) = 0; % Convert to time exposed (make burial periods zero)
      cum_t_expo = fliplr(cumsum(fliplr(t_expo))); % Calculate cumulative exposure time (from present)
  end
  
  % Append reduced time, logical, interval time and cumulative exposure
  % time (add present-day exposure for consistency)
  exposed_or_not.reduced_logical = logical([reduced_logical 1]);
  exposed_or_not.reduced_time = [reduced_time 0];
  exposed_or_not.reduced_values = [reduced_values 1];
  exposed_or_not.reduced_interval_time = [interval_time 1];
  exposed_or_not.cum_t_expo = [cum_t_expo 1];
  
end
    