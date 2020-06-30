function retVal = writeLBData(filename, data, format)
%WRITELBDATA Write Lattice Boltzmann fluid density data to file
%   WRITELBDATA(filename,data,[format]) where data is the matrix to be
%   written to the file named 'filename' in the format specified by
%   'format'. The default format is 'binary'. The format string 'text'
%   outputs the data as plain CSV text. Any other string will write the 
%   data in the binary format. 
%
%   See also readLBData
%
%   Copyright 2009

if ( length(data) == 0 || length(filename) == 0)
    retVal = -1;
    return;
end

if ( nargin < 3 )
    format = 'binary'; % default output format
    if( strcmpi(filename(end-6:end),'.lbData') )
      format = 'text';
    end
end

format = strtrim(format);

if ( strcmpi(format,'text') == 1 )
    dlmwrite(filename,data);
else % assumed to be binary
   retVal = writeLBDataBin(data, filename);
   return;
end
end

function retVal = writeLBDataBin(data, filename)
[fid errmsg] = fopen(filename, 'wb');

if ( fid == -1 )
    display(errmsg);
    retVal = 1;
end

% Find out the number of bytes used by the velocities
sampleVel = data(1);
dataTypeInfo = whos('sampleVel');
floatSize = dataTypeInfo.bytes;

header = sprintf('uint32_t nX*nY*nZ, uint32_t floatSize, float f[19*nX*nY*nZ], float size = %d', floatSize);
headerSize = length(header);
fwrite(fid,headerSize,'uint16');
fwrite(fid,header,'char');

nXnYnZ = size(data,1);
fwrite(fid,nXnYnZ,'uint32');
fwrite(fid,floatSize,'uint32');

fwrite(fid,data,'double');
fclose(fid);
retVal = 0;
end
