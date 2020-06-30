function [models params] = modelQuery(executable);
% MODELQUERY 	Report on supported lattice-Boltzmann models.
% 	[MODELS PARAMS] = MODELQUERY(EXECUTABLE) returns the id strings MODELS 
%   and the corresponding parameter key strings PARAMS for the executable    
%	EXECUTABLE.
%
%   Example:
%       % list the supported material models
%		models = modelQuery('generalLBSimulation');
%       models 
%       % list the supported material models and their parameters
%       [models params] = modelQuery('generalLBSimulation');
%       format long
%       display([models params])
%
%   See also: checkMaterialModel, inputParameterQuery, and 
%     writeMaterialModelFile.
%
%   Copyright 2009

%system('mkfifo lbMatlabPipe'); % create a named pipe 'lbMatlabPipe'
%system('touch lbMatlabPipe'); % create a text file 'lbMatlabPipe' - safer ??


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check for existence of executable
%

x = dir(executable);

% FIXME - This is system specific and a little clunky. 
if(isempty(x))
  % executable not found in current directory - look for it in the system's search path
  runString = ['which ',executable, '> lbMatlabPipe &'];
  system(runString);
  
  fid = fopen('lbMatlabPipe','r');
  tline = fgetl(fid);
  fclose(fid);
  system('rm lbMatlabPipe');
  
  if(isempty(dir(tline)))
    % executable not found in system search path
    disp(['Error, could not find executable "',executable,'".'])
    return
  else   
    % executable is in system search path
    runString = [executable,' --materialModels > lbMatlabPipe &'];
  end
else  
  runString = ['./',executable,' --materialModels > lbMatlabPipe &'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get the models and their parameters
%

% Run the executable and pipe the material models to lbMatlabPipe.
system(runString);

fid = fopen('lbMatlabPipe','r');

models = textscan(fid,'%s %s','Delimiter','{}','MultipleDelimsAsOne',1);

fclose(fid); 

system('rm lbMatlabPipe');

params = strtrim(models{2});
models = strtrim(models{1});


