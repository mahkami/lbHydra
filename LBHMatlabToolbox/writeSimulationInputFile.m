function writeSimulationInputFile(simulationInputFile,simulationInput)
% WRITESIMULATIONINPUTFILE Writes simulation input parameters to a file.
%   writeSimulationInputFile(simulationInputFile,simulationInput) 
% 	    simulationInputFile: The name of the file.
% 	    simulationInput: Array of matlab structures with key fields and string or numerical values, i.e.
%  	 		simulationInput.key =  'value';
%   		simulationInput.key =  'value % comment';
%           simulationInput.key =  value;    % with a numerical value
%
%	See also: checkSimulationInput, inputParameterQuery, and writeMaterialModelFile. 
%
%   Copyright 2009

fid = fopen(simulationInputFile,'w');

keyNames = fieldnames(simulationInput);
for i = 1:length(keyNames)
    keyNameStr = keyNames{i};
    valueStr = getfield(simulationInput,keyNameStr);
    if(~isstr(valueStr))
      valueStr =  num2str(valueStr,9); % convert to string with 9 digits of precision
    end
    fprintf(fid,'%s=%s \n', keyNameStr, valueStr);
end

fclose(fid);
