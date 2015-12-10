function data = load_binListArray(filename, type)
% function data = load_binListArray(filename, type)
%
%  Load an array of lists from a cpplib binary file. Expected file format is number_of_lists, array of ints with list sizes, list elements.
%

fid=fopen(filename,'r');  if fid<=0 fprintf(1,'Error opening file %s!\n',filename); data=[]; return; end;

numLists  = fread(fid,1,'int32');      data = cell(numLists,1);
index = fread(fid,numLists,'int32');   elems = fread(fid,sum(index),type);
elemIdx = 1;
for i=1:numLists
    data{i} = elems(elemIdx:elemIdx+index(i)-1)';
    elemIdx=elemIdx+index(i);
end;

fclose(fid);
