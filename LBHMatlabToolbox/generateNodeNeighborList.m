function [nodeNbrs nodeMaterialModelMap] = generateNodeNeighborList(voxels,c)
%  GENERATENODENEIGHBORLIST		Create a matrix of node neighbors.
%
%    generateNodeNeighborList(VOXELS) uses the 3D array VOXELS to create 
%    a list of node neighbors for the D3Q19 lattice. Elements of VOXELS 
%    less than or equal to 0 are considered holes and are removed from 
%    the list of nodes.
%
%    generateNodeNeighborList(VOXELS,C) uses the lattice directions given 
%    in C in place of the default D3Q19 lattice directions.
%
%    [nodeNbrs,nodeMaterialModelMap] = generateNodeNeighborList(VOXELS)
%    [nodeNbrs,nodeMaterialModelMap] = generateNodeNeighborList(VOXELS,C)
%    creates both a list of the node neighbors and a material model map,
%    which stores the node coordinates, the node and the material model
%    i.e., nodeMaterialModelMap = [X,Y,Z,Node,MaterialModel].
%
%    See also generateNodeMaterialModelMap
%  	 Copyright 2009

%TODO: Extend to support multiple material models
dims = size(voxels);

if ( size(dims) < 3 ) %for compatibility in case of voxels less than 3 dimensions
  if ( size(dims) < 2 )
      if ( size(dims) < 1 ) 
          dims(1) = 1;
      end
      dims(2) = 1;
  end
  dims(3) = 1;
end

if nargin < 2
 % set node connectivity to d3q19 by default
 c = d3q19LatticeDirections();
end

[x,y,z]  = ndgrid([0:(dims(1)-1)],[0:(dims(2)-1)],[0:(dims(3)-1)]);

nodeNbrs = [];

coords = [x(:),y(:),z(:),voxels(:)];

clear x y z;

for i = 1:length(c)

  x = coords(:,1) + c(i,1);
  y = coords(:,2) + c(i,2);
  z = coords(:,3) + c(i,3);
  
  x = mod(x,dims(1))+1;
  y = mod(y,dims(2))+1;
  z = mod(z,dims(3))+1;

  ind = sub2ind(dims,x,y,z);

  nodeNbrs = [nodeNbrs,ind];

end

% remove neighbors according to voxel list
removeList = find(voxels <= 0);

if ( size(removeList,1) > 0 ) % if there are holes
  coords(removeList,:) = [];
  nodeNbrs(removeList,:) = [];

  ind = ismember(nodeNbrs,removeList);
  nodeNbrs(ind) = 0;

  % renumber nodes without gaps
  newNums = cumsum(~~voxels(:)); % we need 0 for holes and 1 for nodes

  ind = find(nodeNbrs(:) > 0);
  nodeNbrs(ind) = newNums(nodeNbrs(ind));
end

% automatic bounceback - set all 0's to row numbers
%[i j] = ind2sub(size(nodeNbrs),ind);
%nodeNbrs(ind) = i;

% generate nodeMaterialModelMap
if nargout > 1   
  nodeNum = cumsum(coords(:,4)>0); % only number nodes that are not holes
  indx = find(coords(:,4)>0);
  nodeMaterialModelMap = [coords(indx,1:3),nodeNum(indx),coords(indx,4)]; % ==> [x y z nodeId materialModel]
end

end
