%
% sample_data = norm_concs_1014(sample_data)
%
% Normalises nuclide concentrations of given samples.
%
% sample_data is a required struct containing necessary sample details, 
% initially created using get_data.m.
%
% Output is the input with the addition of production rates and normalised
% concentrations for each sample and nuclide.
%
%
%%

function sample_data = norm_concs_1014(sample_data)

  % Check inputs
  if (nargin ~= 1)
      error('norm_concs_1014 has wrong number of inputs!');
  end
  if isempty(sample_data.cover.z)
      sample_data.cover.z = 0;
  end
  

  s = sample_data.s;
  
  for a = 1:length(s)
      
      % Compute depths
      this_top_z = s{a}.top_z_gcm2 + sample_data.cover.z;
      this_bottom_z = s{a}.bottom_z_gcm2 + sample_data.cover.z;
      
      
      if s{a}.nuclide10 == 1
          
          % Get sample parameters
          sf = sample_data.sf1026{a};
          cp = sample_data.cp1026{a};
          
          % Define production rate function
          %pr_func = @(z,t) PR_Z(z,sample_data.pp,sf,cp,10);
          pr_func = @(z) PR_Z(z,sample_data.pp,sf,cp,10);
          
          % Integrate over depth for each sample
          Integrated_PR = zeros(1,length(this_top_z));
          for b = 1:length(this_top_z)
              % Integral of production rate in interval divided by interval thickness.
              Integrated_PR(b) = integral(pr_func,this_top_z(b),this_bottom_z(b),'RelTol',1e-3,'AbsTol',1e-3) ./ (this_bottom_z(b) - this_top_z(b));
          end
          
          % Average by mineral weight
          this_weight = s{a}.weight10;
          sample_data.s{a}.PR_10 =  (sum(Integrated_PR.*this_weight)) ./ sum(this_weight);
          
          % Normalise nuclide concentration
          sample_data.s{a}.norm_N10 = s{a}.N10 ./ sample_data.s{a}.PR_10;
          sample_data.s{a}.norm_dN10 = s{a}.dN10 ./ sample_data.s{a}.PR_10;
      
      end
          
      if s{a}.nuclide14 == 1
          
          % Get sample parameters
          sf = sample_data.sf14{a};
          cp = sample_data.cp14{a};
          
          % Define production rate function
          pr_func = @(z,t) PR_Z(z,sample_data.pp,sf,cp,14);
          
          % Integrate over depth for each sample
          Integrated_PR = zeros(1,length(this_top_z));
          for b = 1:length(this_top_z)
              % Integral of production rate in interval divided by interval thickness.
              Integrated_PR(b) = integral(pr_func,this_top_z(b),this_bottom_z(b),'RelTol',1e-3,'AbsTol',1e-3) ./ (this_bottom_z(b) - this_top_z(b));
          end
          
          % Average by mineral weight
          this_weight = s{a}.weight14;
          sample_data.s{a}.PR_14 =  (sum(Integrated_PR.*this_weight)) ./ sum(this_weight);
          
          % Normalise nuclide concentration
          sample_data.s{a}.norm_N14 = s{a}.N14 ./ sample_data.s{a}.PR_14;
          sample_data.s{a}.norm_dN14 = s{a}.dN14 ./ sample_data.s{a}.PR_14;
          
      end
      
  end
  
end
