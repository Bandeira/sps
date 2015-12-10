function [pairs, aligns] = load_resultsPA(filename,separator)
% function [pairs, aligns] = load_resultsPA(filename,separator)

% pairs = int16([]); aligns = single([]);
pairs = uint32([]); aligns = single([]);
fid = fopen(filename,'r'); if fid<=0 fprintf(1,'ERROR opening %s\n',filename); return; end;

numElems = str2double(fgetl(fid));
aligns = single(zeros(numElems,3));
% pairs = int16(zeros(numElems,2));
pairs = uint32(zeros(numElems,2));
for i=1:numElems
    line = fgetl(fid);
    [str,line] = strtok(line,separator);   pairs(i,1) = str2double(str);
    [str,line] = strtok(line,separator);   pairs(i,2) = str2double(str);
    [str,line] = strtok(line,separator);   aligns(i,1) = str2double(str);
    [str,line] = strtok(line,separator);   aligns(i,2) = str2double(str);
    [str,line] = strtok(line,separator);   aligns(i,3) = str2double(str);
    aligns(i,4) = str2double(line(2:size(line,2)));
end;

fclose(fid);
