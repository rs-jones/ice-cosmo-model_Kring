%
% exposed_or_not = time_series_logical(time_series,threshold_frac,model_time,data_use)
% exposed_or_not = time_series_logical(time_series,threshold_frac,model_time,data_use,model_interval,frac_buried)
%
% Creates a logical of time intervals that are either exposed or buried
% based on time-series data and a threshold (1=exposed, 0=buried).
%
% Output is to be used to model nuclide concentrations through time.
%
% time_series should be two columns (time in years before present, and
% corresponding time series values).
%
% threshold_frac should be a fraction of the maximum time_series value 
% (0-1). It is used to determine a threshold value for computing exposure 
% and burial intervals.
%
% data_use should either be 'normal' or 'opposite'. If 'opposite' is set, 
% exposure intervals defined by threshold_frac become burial intervals, 
% and vice versa. Default is normal.
%
% model_interval (in years) is optional. Default is 1000 years.
%
% frac_buried is the fraction of an exposure period that a sample is buried
% (optional).
%
% The interpolated time and values from time_series, threshold, interval 
% time, computed logical, and reduced times and logical are exported to a 
% single structure.
%
%
%%

function exposed_or_not = time_series_logical(time_series,threshold_frac,t_start,data_use,model_interval,frac_buried)
  
  % Check inputs
  if (nargin > 6)
      error('time_series_logical has too many inputs!');
  end
  if (numel(time_series.data(1,:)) ~= 2)
      error('time_series requires two columns (time and corresponding values)!');
  end
  if (time_series.data(1,1) > time_series.data(2,1))
      error('time in time_series needs to be in years before present with the most recent (e.g. 0 yrs) first!');
  end
  if (threshold_frac < 0 && threshold_frac > 1)
      error('threshold_frac must be between 0 and 1!');
  end
  if (strcmp(data_use,'normal') | isempty(data_use))
      opposite = 0;
  elseif (strcmp(data_use,'opposite'))
      opposite = 1;
  else
      error('data_use should be "normal" or "opposite"!');
  end
  if (nargin < 5) || isempty(model_interval)
      model_interval = 1000;
  end
  if (model_interval ~= 1000 && model_interval ~= 100 && model_interval ~= 10)
      error('model_interval must be 10, 100 or 1000 years!');
  end
  if (nargin < 6) || isempty(frac_buried)
      frac_buried = 0;
  end
  
  
  if model_interval == 1000
      r_fac = -3;
  elseif model_interval == 100
      r_fac = -2;
  elseif model_interval == 10
      r_fac = -1;
  end
  
  t_start_a = round(t_start,r_fac);
  
  
  % Add time series name and units to output
  exposed_or_not.name = time_series.name;
  exposed_or_not.short_name = time_series.short_name;
  exposed_or_not.unit = time_series.unit;
  
  % Convert time_series to model_interval spacing
  time = 0:model_interval:t_start_a;  % Create new time vector
  ts_values = interp1(time_series.data(:,1),time_series.data(:,2),time);
  
  % For any time intervals where the value is NaN, make value the same as
  % the previous time interval
  ts_values = fliplr(ts_values); % Flip to work forward in time
  for i = 1:(numel(ts_values)-1) % Avoid oldest time interval
      if isnan(ts_values(i+1))
          ts_values(i+1) = ts_values(i);
      end
  end
  ts_values = fliplr(ts_values); % Flip back
  
    
  % Compute threshold value
  threshold = min(time_series.data(:,2)) + ((max(time_series.data(:,2))-min(time_series.data(:,2))) * threshold_frac);
  
  % Determine intervals of exposure and burial using threshold and create logical
  % default: >threshold = exposure, <threshold = burial
  if opposite == 1
      exposed_or_not.logical = logical(ts_values < threshold); 
  else
      exposed_or_not.logical = logical(ts_values > threshold);
  end
  
  % Flip data so that oldest time is first
  time = fliplr(time);
  ts_values = fliplr(ts_values);
  exposed_or_not.logical = fliplr(exposed_or_not.logical);
  
  % Append time, values, threshold and interval
  exposed_or_not.time = time;
  exposed_or_not.values = ts_values;
  exposed_or_not.interval_time = zeros(1,numel(exposed_or_not.logical)) + model_interval;
  exposed_or_not.threshold = threshold;
  
  
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
      end
      t_expo = interval_time; t_expo(~reduced_logical) = 0; % Convert to time exposed (make burial periods zero)
      cum_t_expo = fliplr(cumsum(fliplr(t_expo))); % Calculate cumulative exposure time (from present)
  end
  
  % Append reduced time, logical, interval time and cumulative exposure 
  % time (add present-day exposure for consistency)
  exposed_or_not.reduced_logical = logical([reduced_logical 1]);
  exposed_or_not.reduced_time = [reduced_time 0];
  exposed_or_not.reduced_interval_time = [interval_time 1];
  exposed_or_not.cum_t_expo = [cum_t_expo 1];
  
  % Get time of initial exposure
  initial_idx = find(exposed_or_not.reduced_logical,1);
  time_initial = exposed_or_not.reduced_time(initial_idx);
  if ~isempty(time_initial)
      exposed_or_not.time_initial = time_initial;
  else
      exposed_or_not.time_initial = 0;
  end
  
end
    