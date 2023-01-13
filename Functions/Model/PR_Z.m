%
% Ptotal = PR_Z(z_in,pp,sf,cp,nuclide)
%
% Calculates total production rates (from spallation and muon pathways) for
% given depths based on computed sample parameters, using CronusCalc 
% functions.
%
% z_in should be a depth or array of depths (in g/cm^2).
%
% pp is a struct of physical constant, generated using get_pars.m.
%
% sf is a struct of scaling factors, generated for each sample using 
% get_pars.m.
%
% cp is a struct of computed parameters, generated for each sample using 
% get_pars.m.
%
% nuclide should 10, 26 or 14, corresponding to 10Be, 26Al, or 14C. 
%
% Output is the total production rate at given depths for the given 
% nuclide, for all samples.
%
%
%%

function Ptotal = PR_Z(z_in,pp,sf,cp,nuclide)

  % Check inputs
  if (nargin ~= 5)
      error('PR_Z has wrong number of inputs!');
  end

  if isvector(z_in)
      z = z_in;
  else
      z = reshape(z_in,1,numel(z_in)); % Unwrap matrix input
  end

  if nuclide == 10 || nuclide == 26
  
      [PtotalBe,PtotalAl,~,~,~,~] = prodz1026(z,pp,sf,cp);
      
      if nuclide == 10
          Ptotal = PtotalBe;
      elseif nuclide == 26
          Ptotal = PtotalAl;
      end
      
  elseif nuclide == 14
      [Ptotal,~,~] = prodz14(z,pp,sf,cp);

  end
  
  if ~isvector(z_in)
      Ptotal = reshape(Ptotal,size(z_in,1),size(z_in,2)); % Re-wrap matrix input
  end

end
