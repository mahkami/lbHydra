function F = calculateLBData(densities,velocities,weights,directions);
% CALCULATELBDATA	Create lattice-Boltzmann density distributions.  
% 	 CALCULATELBDATA(RHO) returns lattice-Boltzmann density data 
%    calculated from the vector of fluid densities "RHO".

%	CALCULATELBDATA(RHO,V) returns lattice-Boltzmann densities based 
%   on from the fluid densities "RHO" and velocities "V".
%
%   CALCULATELBDATA(RHO,[],WEIGHTS) uses WEIGHTS in place of the 
%   default D3Q19 lattice weights. 
%
%   CALCULATELBDATA(RHO,V,WEIGHTS,DIRECTIONS) uses WEIGHTS and 
%   DIRECTIONS in place of the default D3Q19 lattice directions and 
%   weights. 
%
%   See also calculateDensities, d2q9LatticeWeights,
%   d2q9LatticeDirections, d3q19LatticeWeights, 
%   d3q19LatticeDirections
%
%	Copyright 2010

% Assume d3q19 by default
if nargin < 4
  directions =  d3q19LatticeDirections;
end
if nargin < 3
  weights = d3q19LatticeWeights;
end
if nargin < 2
  velocities = [];
end

F = densities(:)*weights(:)';
numDirs = length(weights);

if( ~isempty(velocities) )
  threeXDotProds = 3*velocities*directions';

  vv = sum(velocities.^2,2);

  F = F.*(1.0 + threeXDotProds.*(1.0 + threeXDotProds/2.0) -  1.5*vv(:,ones(numDirs,1)) );
end
