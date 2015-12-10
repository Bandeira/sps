function data = sn_load_binarray(filename, type)
% function data = sn_load_binarray(filename, type)

fid=fopen(filename,'r');  if fid<=0 fprintf(1,'Error opening file %s!\n',filename); data=[]; return; end;

dim  = fread(fid,2,'int32');
data = fread(fid,dim(1)*dim(2),type);
data = reshape(data,dim(2),dim(1))';

fclose(fid);