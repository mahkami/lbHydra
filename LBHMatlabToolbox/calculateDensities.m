function rho = calculateDensities(A)
% CALCULATEDENSITIES		Calculate densities from lattice-Boltzmann data.
%	CALCULATEDENSITIES(LBDATA) Returns the density of the lattice-Boltzmann data 
%   in LBDATA.
%
%   See also calculateVelocities.
%	Copyright 2009
rho = sum(A,2);
