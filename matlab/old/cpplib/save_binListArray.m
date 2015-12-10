function save_binListArray(filename, data, dataType)
%
%  Saves a vector of lists to a file.
%
%  data - cell vector where data{i} is a vector of dataType elements.
%

fid=fopen(filename,'w');  if fid<=0 fprintf(1,'ERROR opening %s!\n',filename); return; end;
numLists=size(data,1);   index=zeros(numLists,1);
for i=1:numLists index(i)=length(data{i}); end;
fwrite(fid,numLists,'uint');
fwrite(fid,index,'uint');
for i=1:numLists fwrite(fid,data{i},dataType); end;
fclose(fid);
