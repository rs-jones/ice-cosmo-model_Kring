%
% sample_data = get_data_1014(input_name)
% 
% Loads input .xlsx, extracts sample names and concentrations, sorts 
% data for CRONUScalc functions, and exports to single output.
%
% This version also sorts data for plotting multiple nuclide 
% concentrations.
%
% For CRONUScalc calculations, the erosion rate and inheritance is set to 
% zero. Attenuation lengths are automatically determined using the Sato 
% model, based on the location information of each sample.
%
% Input data must include the following information (Note difference to 
% required data input for get_data.m):
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
% Written by Richard Selwyn Jones, Monash University
% richard.s.jones@monash.edu
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function sample_data = get_data_1014(input_name)
  
  % Load input data
  [in_data,in_txt,~] = xlsread(strcat(input_name,'.xlsx'));
  
  % Check inputs
  if (numel(in_txt(1,:)) > 25)
      error('sample input data has too many fields!');
  end
  if (numel(in_txt(1,:)) > 18 && numel(in_txt(1,:)) < 25)
      error('sample input data seems to have exposure ages but fields are missing!');
  end
  
  if (numel(in_txt(1,:)) > 18)
      ages = in_data(:,18:23);
      scaling = in_txt(2:end,end);
  else
      ages = [];
      scaling = [];
  end
  
  
  % Get number of samples (and samples with 10-Be and 14-C concs.)
  n_samples = numel(in_data(:,1));
  
  % Get sample names
  names = in_txt(2:end,1)';
  
  
  % If concentrations are NaN, make zero
  for a = 1:n_samples
      if (isnan(in_data(a,9)))
          in_data(a,9) = 0;
          in_data(a,10) = 0;
      end
      if (isnan(in_data(a,11)))
          in_data(a,11) = 0;
          in_data(a,12) = 0;
      end
  end

  % Determine thickness or sample depths if not all set
  for b = 1:n_samples
      if (isempty(in_data(b,6)) || in_data(b,6) == 0 || isnan(in_data(b,6)))
          in_data(b,6) = in_data(b,14) - in_data(b,13);
      end
      if (isempty(in_data(b,14)) || isnan(in_data(b,14)))
          in_data(b,14) = in_data(b,15) - in_data(b,6);
      end
      if (isempty(in_data(b,15)) || isnan(in_data(b,15)))
          in_data(b,15) = in_data(b,14) + in_data(b,6);
      end
  end
  
 
  % Sort details for each sample and nuclide (currently only 10Be and 14C)
  logical_10 = any(in_data(:,9),2)';
  logical_14 = any(in_data(:,11),2)';
    
  NN = [any(logical_10) any(logical_14)];
  
  % Determine if sample details should be combined for each nuclide measurement, and sort if necessary
  if NN(1)
      N = in_data(logical_10,9:10);
      in_data_N = in_data(logical_10,:);
      sorted_10 = sort_data_1014(in_data_N,N,names);
  else
      sorted_10 = [];
  end
  if NN(2)
      N = in_data(logical_14,11:12);
      in_data_N = in_data(logical_14,:);
      sorted_14 = sort_data_1014(in_data_N,N,names);
  else
      sorted_14 = [];
  end
  
  logical_1014 = any(logical_10' & logical_14'); % Determine whether samples have both nuclide measurements
  
  % Get data for samples with multiple nuclide measurements
  if logical_1014
      sorted_names = sorted_10.names;
      sorted_data = sorted_10.data;
      duplicate_logical = sorted_10.duplicates;
      duplicates_10 = sorted_10.duplicates;
      duplicates_14 = sorted_14.duplicates;
  % Get data for samples with a single nuclide
  elseif NN(1) & ~NN(2)
      sorted_names = sorted_10.names;
      sorted_data = sorted_10.data;
      duplicate_logical = sorted_10.duplicates;
      duplicates_10 = sorted_10.duplicates;
  % Combine sorted data for nuclides if samples do not have multiple nuclide measurements
  else %~logical_1026
      sorted_names = [sorted_10.names; sorted_14.names];
      sorted_data = [sorted_10.data; sorted_14.data];
      duplicate_logical = [sorted_10.duplicates sorted_14.duplicates];
      duplicates_10 = [sorted_10.duplicates zeros(1,length(sorted_14.duplicates))];
      duplicates_14 = [zeros(1,length(sorted_10.duplicates)) sorted_14.duplicates];
  end
  
  % Calculate number of sorted samples
  n_sorted = numel(sorted_data(:,1));
  logical_sorted_10 = any(sorted_data(:,9),2)';
  logical_sorted_14 = any(sorted_data(:,11),2)';
  
  
  % Collate sample details
  for c = 1:n_sorted
      
      sample_data.s{c}.name = sorted_names(c);
      
      if logical_sorted_10(c)
          sample_data.s{c}.nuclide10 = 1;
          sample_data.s{c}.N10 = sorted_data(c,9);
          sample_data.s{c}.dN10 = sorted_data(c,10);
          sample_data.s{c}.weight10 = sorted_data(c,15);
      else
          sample_data.s{c}.nuclide10 = 0;
      end
      if logical_sorted_14(c)
          sample_data.s{c}.nuclide14 = 1;
          sample_data.s{c}.N14 = sorted_data(c,11);
          sample_data.s{c}.dN14 = sorted_data(c,12);
          sample_data.s{c}.weight14 = sorted_data(c,16);
      else
          sample_data.s{c}.nuclide14 = 0;
      end
      
      sample_data.s{c}.top_z = sorted_data(c,13);
      sample_data.s{c}.bottom_z = sorted_data(c,14);
      
  end
  
  % Add all details for each nuclide measurement (if samples were combined)
  if any(duplicate_logical)
     if any(duplicates_10)
         dup_s = find(duplicates_10);
         for d = 1:length(dup_s)
             s = dup_s(d);
             sample_data.s{s}.name = sorted_10.dup_names{d};
             sample_data.s{s}.top_z = sorted_10.dup_top_z{d};
             sample_data.s{s}.bottom_z = sorted_10.dup_bottom_z{d};
             sample_data.s{s}.weight10 = sorted_10.dup_weight10{d};
             sample_data.s{s}.weight14 = sorted_10.dup_weight14{d};
         end
     end
     if any(duplicates_14)
         dup_s = find(duplicates_14);
         for d = 1:length(dup_s)
             s = dup_s(d);
             sample_data.s{s}.name = sorted_10.dup_names{d};
             sample_data.s{s}.top_z = sorted_14.dup_top_z{d};
             sample_data.s{s}.bottom_z = sorted_14.dup_bottom_z{d};
             sample_data.s{s}.weight10 = sorted_10.dup_weight10{d};
             sample_data.s{s}.weight14 = sorted_10.dup_weight14{d};
         end          
     end
  end
  
  
  % Calculate depths in g/cm^2 for each sample
  for e = 1:n_sorted
      
      sample_data.s{e}.top_z_gcm2 = sample_data.s{e}.top_z .* sorted_data(e,7);
      sample_data.s{e}.bottom_z_gcm2 = sample_data.s{e}.bottom_z .* sorted_data(e,7);
      
      % Combine top depths in g/cm^2
      top_depth_gcm2(e) = mean(sample_data.s{e}.top_z_gcm2);
  
  end
  
  
  % Assume certain sample details are unknown or zero
  e_rate = zeros(n_sorted,1); % Erosion rate (mm/kyr)
  inh10 = zeros(n_sorted,1);  % 10Be inheritance (atoms/g)
  inh14 = zeros(n_sorted,1);  % 14C inheritance (atoms/g)
%   for i = 1:n_sorted
%       min14(i,1) = 'q'; % 14C mineral phase ('q'=quartz)
%   end
  
  % Initially set attenuation length to zero
  init_L = zeros(n_sorted,1);
  
  % Combine necessary data for CronusCalc
  sample_data.CC.Be10 = [sorted_data(:,[1:4,6:8]),e_rate,sorted_data(:,9),zeros(size(sorted_data(:,9))),inh10,inh14,init_L,top_depth_gcm2',sorted_data(:,17)];
  sample_data.CC.C14 = [sorted_data(:,[1:4,6:8]),e_rate,sorted_data(:,9),inh14,init_L,top_depth_gcm2',sorted_data(:,17)];
  
  % Determine atmospheric pressure (if unknown)
  sample_data.CC.Be10 = atm_pressure(sorted_data(:,1),sorted_data(:,2),sample_data.CC.Be10);
  sample_data.CC.C14 = atm_pressure(sorted_data(:,1),sorted_data(:,2),sample_data.CC.C14);

  % Determine attenuation lengths (g/cm^2) based on Sato model
  for d = 1:length(sample_data.CC(:,1))
      L(d,1) = attenuationlength(sample_data.CC.Be10(d,1),sample_data.CC.Be10(d,2),sample_data.CC.Be10(d,3),sample_data.CC.Be10(d,4));
  end
  sample_data.CC.Be10(:,13) = L; % Add attenuation lengths
  sample_data.CC.C14(:,11) = L; % Add attenuation lengths
  
  
  % Create uncertainties struct
  
  % Create a standard elevation uncertainty if not provided
  for eu = 1:length(sorted_data(:,3))
      if sorted_data(eu,3) == 0 || isnan(sorted_data(eu,3))
          elev_uncert(eu) = 5; % 5 metres
      else
          elev_uncert(eu) = sorted_data(eu,3);
      end
  end
  uncert_10 = zeros(size(sample_data.CC.Be10)); uncert_14 = uncert_10; % Assume zero uncertainty for now
  uncert_10(:,9) = sorted_data(:,10); % Add nuclide concentration uncertainties
  uncert_14(:,9) = sorted_data(:,12); % Add nuclide concentration uncertainties
  sample_data.CC_uncert.Be10 = uncert_10;
  sample_data.CC_uncert.C14 = uncert_14;
  
  % Create ages struct
  if ~isempty(ages)
      if NN(1)
          sample_data.ages.Be10 = ages(:,1:3);
      end
      if NN(2)
          sample_data.ages.C14 = ages(:,4:6);
      end
      sample_data.ages.scaling_model = scaling{1}; % Assume the scaling used is the same for all samples
  end
  
  % Export logical of nuclides
  sample_data.logical_10 = logical_sorted_10;
  sample_data.logical_14 = logical_sorted_14;
  
  % Export sample position
  sample_data.position = sorted_data(:,5)';
  
  % Create default cover depth
  sample_data.cover.z = 0;
  
  
end
