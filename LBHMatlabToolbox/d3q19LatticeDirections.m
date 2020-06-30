function c = d3q19LatticeDirections();
% D3Q19LATTICEDIRECTIONS	Neighbor node directions for the D3Q19 lattice-Boltzmann lattice.
%   D3Q19LATTICEWEIGHTS() returns a 19x3 matrix of vectors for the 
%   D3Q19 lattice-Boltzmann lattice. 
%
%   See also d3q19LatticeWeights
%	Copyright 2009


c = [  0,0 ,0
      -1, 0, 0 
      1, 0, 0
      0,-1, 0
      0, 1, 0
      0, 0, -1
      0, 0, 1
     -1,-1, 0
      1,-1, 0
     -1, 1, 0
      1, 1, 0 
     -1,0, -1
      1,0,-1
     -1,0, 1
      1,0, 1 
      0,-1,-1
      0, 1,-1
      0,-1, 1
      0, 1, 1 ];

