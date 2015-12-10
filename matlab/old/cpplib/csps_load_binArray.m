function data = load_binArray(filename, type)
% function data = load_binArray(filename, type)
%
%  Load a binary array from cpplib. Expected file format is number_of_lines,number_of_cols,data.
%

fid=fopen(filename,'r');  if fid<=0 fprintf(1,'Error opening file %s!\n',filename); data=[]; return; end;

dim  = fread(fid,2,'int32');
data = fread(fid,dim(1)*dim(2),type);
data = reshape(data,dim(2),dim(1))';

fclose(fid);