function c = d3q7LatticeDirections();
% D3Q7LATTICEDIRECTIONS	Neighbor node directions for the D3Q7 lattice-Boltzmann lattice.
%   D3Q7LATTICEWEIGHTS() returns a 7x3 matrix of vectors for the 
%   D3Q7 lattice-Boltzmann lattice. 
%
%   See also d3q7LatticeWeights
%	Copyright 2009


c = [  0,0 ,0
      -1, 0, 0 
      1, 0, 0
      0,-1, 0
      0, 1, 0
      0, 0, -1
      0, 0, 1];

