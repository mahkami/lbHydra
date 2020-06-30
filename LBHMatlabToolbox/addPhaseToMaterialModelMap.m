function newMaterialModelMap = addPhaseToMaterialModelMap(nodeMaterialModelMap,voxels)
%  ADDPHASETOMATERIALMODELMAP	Add a phase to a node-material-model-map
%    addPhaseToMaterialModelMap(MM_MAP,VOXELS) Adds the material model ids 
%    in VOXELS to the material model map MM_MAP as a new phase and returns 
%    the result.
%
%    See also: generateNodeNeighborList
%
%    Copyright 2009

% flatten voxels;
voxels =voxels(:);

% remove neighbors according to voxel list
removeList = find(voxels <= 0);
voxels(removeList) = [];

if(length(voxels) ~= length(nodeMaterialModelMap))
  display('Error: number of material models in new phase not equal to number of nodes in material model map.')
  return
end

% append new material models to existing material model map
newMaterialModelMap = [nodeMaterialModelMap,voxels];

end
