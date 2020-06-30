function fprintMaterialModel(fid,materialModel,mmNumber)
% FPRINTMATERIALMODEL 	Write material model parameters to a file id.
%   fprintMaterialModel(fid,materialModel,mmNumber)
% 		materialModel: is a struct with fields "number", "name" and "params".
%   		materialModels(i).name: lists the name of the ith material model
%   		materialModels(i).params: struct with string fields for each
%           material model parameter, eg.
%           	materialModels(i).params.normalized_viscosity = '0.1667'
%   
%	See also: writeMaterialModelFile, writeSimulationInputFile
%
%   Copyright 2009

  if(nargin > 2) 
    fprintf(fid,'number=%d,name=%s,params=',mmNumber,materialModel.name);
  elseif(isfield( materialModel,'number'))
    fprintf(fid,'number=%d,name=%s,params=',materialModel.number,materialModel.name);
  else
    % warning: number not defined
    display('Warning material model number not defined, set to 0')
    fprintf(fid,'number=0,name=%s,params=',materialModel.name);
  end
  
  if( isfield( materialModel,'params') & length(materialModel.params)>0)
    paramNames = fieldnames(materialModel.params);
    for j = 1:length(paramNames)
      paramNameStr = paramNames{j};
      paramStr = getfield(materialModel.params,paramNameStr);
      if(~isstr(paramStr))
        paramStr = sprintf('%.9g ',paramStr); % convert to string with 9 digits of precision
      end
      
      if (j > 1) 
        fprintf(fid,';');
      end
      fprintf(fid,'%s:%s',paramNameStr,paramStr);  
    end
  else 
    fprintf(fid,'none');
  end
  fprintf(fid,'\n'); 
