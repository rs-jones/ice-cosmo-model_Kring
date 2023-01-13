%
% N_data = N_to_plot_1014(sample_data)
% N_data = N_to_plot_1014(sample_data,predN)
% 
% Sorts measured (sample_data) or predicted (predN) nuclide concentrations 
% for plotting. If predN is excluded then concentrations are sorted for 
% measured data only, if predN is included then only the predicted 
% concentrations are normalised and sorted.
%
% sample_data is a required struct, created using get_data.m. It must
% include the normalised concentrations, calculated using norm_concs.m.
%
% predN is data struct of predicted nuclide concentrations. The summed 
% concentrations derived from the complex history model are used if
% present. Otherwise it uses whatever the predN input is, but must be a
% struct containing at least nuclide variables (e.g. N10).
%
% In cases where sample material has multiple concentration measurements
% for one nuclide and only a single measurement for another nuclide (e.g. 
% 2x 10Be measurements but 1x 14C measurement from 2x core segments), the 
% multiple concentrations (i.e. for 10Be in above example) are mixed. This 
% is determined from which samples have common depths.
%
% Output is used with plot_concs_1014.m. It is a struct of normalised 
% nuclide concentrations for all samples, for each nuclide.
% 
%
%%

function N_data = N_to_plot_1014(sample_data,predN)

  % Check inputs
  if (nargin > 2)
      error('N_to_plot_1014 has wrong number of inputs!');
  end
  
  try
      pred = predN.sum; % Use summed predicted nuclide concentrations
  catch
      if (nargin < 2)
          pred = [];
      else
          pred = predN; % If summed concentrations are missing, use whatever the input is
      end
  end

  
  % Isolate sample data for each nuclide
  if any(sample_data.logical_10)
      data_N10 = sample_data.s(sample_data.logical_10);
      % Get top and bottom depths for each sample measurement
      for a = 1:length(data_N10)
          top_z_10(a) = min(data_N10{a}.top_z);
          bottom_z_10(a) = max(data_N10{a}.bottom_z);
      end
  end
  if any(sample_data.logical_14)
      data_N14 = sample_data.s(sample_data.logical_14);
      % Get top and bottom depths for each sample measurement
      for a = 1:length(data_N14)
          top_z_14(a) = min(data_N14{a}.top_z);
          bottom_z_14(a) = max(data_N14{a}.bottom_z);
      end
  end
  
  
  % Find common depths between nuclides
  if any(sample_data.logical_10) && any(sample_data.logical_14)
      common_top_z_1014 = intersect(top_z_10,top_z_14);
      common_bottom_z_1014 = intersect(bottom_z_10,bottom_z_14);
  end
  
  
  % Find samples to mix for each nuclide
  if any(sample_data.logical_10)
      % Create logical matrix (common depths x sample measurements)
      logical_mix_10 = true(length(common_top_z_1014),length(top_z_10));
      % Get sample depths
      for b = 1:length(common_top_z_1014)
          for c = 1:length(top_z_10)
              % Find if sample depth is within range of common depths, for each sample measurement and common depth;
              % Save to matrix
              logical_mix_10(b,c) = (top_z_10(c) >= common_top_z_1014(b) & top_z_10(c) <= common_bottom_z_1014(b));
          end
      end
  end
  if any(sample_data.logical_14)
      % Create logical matrix (common depths x sample measurements)
      logical_mix_14 = true(length(common_top_z_1014),length(top_z_14));
      % Get sample depths
      for b = 1:length(common_top_z_1014)
          for c = 1:length(top_z_14)
              % Find if sample depth is within range of common depths, for each sample measurement and common depth;
              % Save to matrix
              logical_mix_14(b,c) = (top_z_14(c) >= common_top_z_1014(b) & top_z_14(c) <= common_bottom_z_1014(b));
          end
      end
  end
  
  
  % Sort data for plotting
  
  if isempty(pred) % Do for sample data
      
      if any(~all(logical_mix_10,2)) || any(~all(logical_mix_14,2)) % Mix data
          if any(sample_data.logical_10)
              for d = 1:length(common_top_z_1014)
                  log_10 = logical_mix_10(d,:);
                  for e = 1:length(data_N10)
                      norm_N10_all(e) = data_N10{e}.norm_N10;
                      norm_dN10_all(e) = data_N10{e}.norm_dN10;
                      sum_weights_10_all(e) = sum(data_N10{e}.weight10);
                  end
                  norm_N10 = norm_N10_all(log_10);
                  norm_dN10 = norm_dN10_all(log_10);
                  sum_weights_10 = sum_weights_10_all(log_10);
                  % Combine concentrations
                  for f = 1:length(norm_N10)
                      norm_N10_weight(f) = norm_N10(f) .* sum_weights_10(f);
                      norm_dN10_weight_sq(f) = ((norm_dN10(f) .* sum_weights_10(f)) ./ sum(sum_weights_10))^2;
                  end
                  mixed_norm_N10 = sum(norm_N10_weight) ./ sum(sum_weights_10);
                  mixed_norm_dN10 = sqrt(sum(norm_dN10_weight_sq));
                  % Export mixed concentrations
                  N_data.norm_N10(d) = mixed_norm_N10;
                  N_data.norm_dN10(d) = mixed_norm_dN10;
              end
          end
          if any(sample_data.logical_14)
              for d = 1:length(common_top_z_1014)
                  log_14 = logical_mix_14(d,:);
                  for e = 1:length(data_N14)
                      norm_N14_all(e) = data_N14{e}.norm_N14;
                      norm_dN14_all(e) = data_N14{e}.norm_dN14;
                      sum_weights_14_all(e) = sum(data_N14{e}.weight14);
                  end
                  norm_N14 = norm_N14_all(log_14);
                  norm_dN14 = norm_dN14_all(log_14);
                  sum_weights_14 = sum_weights_14_all(log_14);
                  % Combine concentrations
                  for f = 1:length(norm_N14)
                      norm_N14_weight(f) = norm_N14(f) .* sum_weights_14(f);
                      norm_dN14_weight_sq(f) = ((norm_dN14(f) .* sum_weights_14(f)) ./ sum(sum_weights_14))^2;
                  end
                  mixed_norm_N14 = sum(norm_N14_weight) ./ sum(sum_weights_14);
                  mixed_norm_dN14 = sqrt(sum(norm_dN14_weight_sq));
                  % Export mixed concentrations
                  N_data.norm_N14(d) = mixed_norm_N14;
                  N_data.norm_dN14(d) = mixed_norm_dN14;
              end
          end
          
      else % Don't mix
          if any(sample_data.logical_10) && any(sample_data.logical_14)
              for d = 1:length(sample_data.s)
                  N_data.norm_N10(d) = sample_data.s{d}.norm_N10;
                  N_data.norm_dN10(d) = sample_data.s{d}.norm_dN10;
                  N_data.norm_N14(d) = sample_data.s{d}.norm_N14;
                  N_data.norm_dN14(d) = sample_data.s{d}.norm_dN14;
              end
          end
      end
  
      
  else % Do for predicted nuclide concentrations
      
      if any(~all(logical_mix_10,2)) || any(~all(logical_mix_14,2)) % Mix data
          if any(sample_data.logical_10)
              pred_N10 = pred.N10(sample_data.logical_10);
              for d = 1:length(pred_N10)
                  % Normalise nuclide concentrations
                  norm_N10_all(d) = pred_N10(d) ./ data_N10{d}.PR_10;
              end
              for e = 1:length(common_top_z_1014)
                  log_10 = logical_mix_10(e,:);
                  for f = 1:length(data_N10)
                      sum_weights_10_all(f) = sum(data_N10{f}.weight10);
                  end
                  norm_N10 = norm_N10_all(log_10);
                  sum_weights_10 = sum_weights_10_all(log_10);
                  % Combine concentrations
                  for h = 1:length(norm_N10)
                      norm_N10_weight(h) = norm_N10(h) .* sum_weights_10(h);
                  end
                  mixed_norm_N10 = sum(norm_N10_weight) ./ sum(sum_weights_10);
                  % Export mixed concentrations
                  N_data.norm_N10(e) = mixed_norm_N10;
              end
          end
          if any(sample_data.logical_14)
              pred_N14 = pred.N14(sample_data.logical_14);
              for d = 1:length(pred_N14)
                  % Normalise nuclide concentrations
                  norm_N14_all(d) = pred_N14(d) ./ data_N14{d}.PR_14;
              end
              for e = 1:length(common_top_z_1014)
                  log_14 = logical_mix_14(e,:);
                  for g = 1:length(data_N14)
                      sum_weights_14_all(g) = sum(data_N14{g}.weight14);
                  end
                  norm_N14 = norm_N14_all(log_14);
                  sum_weights_14 = sum_weights_14_all(log_14);
                  % Combine concentrations
                  for i = 1:length(norm_N14)
                      norm_N14_weight(i) = norm_N14(i) .* sum_weights_14(i);
                  end
                  mixed_norm_N14 = sum(norm_N14_weight) ./ sum(sum_weights_14);
                  % Export mixed concentrations
                  N_data.norm_N14(e) = mixed_norm_N14;
              end
          end
          
      else % Don't mix
          % Normalise nuclide concentrations
          if any(sample_data.logical_10)
              for d = 1:length(pred.N10)
                  N_data.norm_N10(d) = pred.N10(d) ./ sample_data.s{d}.PR_10;
              end
          end
          if any(sample_data.logical_14)
              for d = 1:length(pred.N14)
                  N_data.norm_N14(d) = pred.N14(d) ./ sample_data.s{d}.PR_14;
              end
          end
      end
  end
    
end
