function data = load_floatArray(filename,numCols)
% function data = load_floatArray(filename,numCols)
%
%  Loads an array of floats with numCols cols from a binary file. 
%  Expected file format:
%    Number of entries - 1 int (4 bytes)
%    numCols * numEntries floats
%

fid = fopen(filename,'r'); if fid<=0 fprintf(1,'Error opening %s!\n',filename); data=[]; return; end;

numEntries = fread(fid,1,'int32'); 
data = fread(fid,numEntries*numCols,'float32')';
data = reshape(data,numCols,numEntries)';

fclose(fid);