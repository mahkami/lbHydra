function writeMaterialModelFile(matModelFile,materialModels)
% WRITEMATERIALMODELFILE	writes material model parameters to a file.
% 	writeMaterialModelFile(matModelFile,materialModels)
% 		matModelFile: The output file string.
%		materialModels: An array of matlab structures with fields "number", 
%       "name" and "params".
%   		materialModels(i).name: Lists the name of the ith material model
%   		materialModels(i).params: Struct with string fields for each 
%           material model parameter, eg.
%           materialModels(i).params.normalized_viscosity = '0.1667'
%
%  See also checkMaterialModel, fprintMaterialModel, materialModelQuery
%           and writeSimulationInputFile.
% 
%  Copyright 2009  

fid = fopen(matModelFile,'w');
for i = 1:length(materialModels)

  if isfield( materialModels(i),'number')
    fprintMaterialModel(fid,materialModels(i));
  else
    fprintMaterialModel(fid,materialModels(i),i);
  end

end
fclose(fid);
