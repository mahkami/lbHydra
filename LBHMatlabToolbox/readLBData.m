function f = readLBData(filename, format)
% READLBDATA Read Lattice Boltzmann fluid density data from file.
%   LBDATA = READLBDATA(FILENAME) returns the lattice boltzmann
%   data in the file FILENAME. The default format is 'binary' 
%   unless FILENAME ends in '.lbData'. 
%       The data is returned as a matrix, with rows equal to the
%   number of nodes and a number of columns equal to the number of 
%   lattice-Boltzmann lattice directions.   
%
%   LBDATA = READLBDATA(FILENAME,FORMAT) reads the data from a 
%   file in the designated format. Supported format strings are
%   'binary' and 'text'. Unrecognized strings will lead to reading 
%   the file in plain text CSV format.
%
%   
%   See also writeLBData.
%
%   Copyright 2009

  if ( nargin < 2 )
    format = 'binary'; % default input format
    if( strcmpi(filename(end-6:end),'.lbData') )
      format = 'text';
    end
  end

  if ( strcmpi(format,'binary') )
    f = readLBBinData(filename);
  else % assumed to be text
    f = dlmread(filename);
  end

end

function f = readLBBinData(filename)
% READLBBINDATA read lattice boltzmann data from a binary file.

f = [];

[fid errmsg] = fopen(filename,'r');

if ( fid == -1 )
    display(errmsg);
    return;
end

n = fread(fid,1,'uint16'); % number of string characters
text = fread(fid,n,'uint8=>char'); % text of string
%display(text');

numNodes = fread(fid,1,'uint32'); % number of nodes

floatSize = fread(fid,1,'uint32')*8; % size of floats

ff = fread(fid,numNodes,['float',int2str(floatSize)]);
while ~feof(fid)
   f = [f,ff];
   ff = fread(fid,numNodes,['float',int2str(floatSize)]);
end

if size(ff,1) == numNodes
  f = [f;ff];
end

fclose(fid);

end
