function data = sn_load_aligns(filename,numCols)
% function data = sn_load_aligns(filename,numCols)

fid = fopen(filename,'r'); if fid<=0 fprintf(1,'Error opening %s!\n',filename); data=[]; return; end;

numEntries = fread(fid,1,'int32'); 
data = fread(fid,numEntries*numCols,'float32')';
data = reshape(data,numCols,numEntries)';
data(:,1:2) = round(data(:,1:2));

fclose(fid);
