function v = calculateVelocities(A,c)
%  CALCULATEVELOCITIES	Calculate velocities from lattice-Boltzmann data.
%    CALCULATEVELOCITIES(LBDATA) calculates the velocities from D3Q19  
%    lattice-Boltzmann data. 
%
%    CALCULATEVELOCITIES(LBDATA,C) uses the lattice directions given 
%    in C in place of the default D3Q19 lattice directions. 
%
%    See also calculateDensities, d3q19LatticeDirections
%	 Copyright 2009

if nargin < 2
  c = d3q19LatticeDirections();
end

rho = sum(A,2);
ind = find(rho<=0);
rho(ind) = 1;

v = A*c;

v = v./rho(:,[1,1,1]);

v(ind,:) = 0;
