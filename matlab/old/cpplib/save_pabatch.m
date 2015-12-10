function save_pabatch(filename, pairs, shiftLists)
%
%  Save job info file of pairwise alignments to be executed by the c++ tools
%
%  File format is
%    Line 1:   NumSpecPairs
%    Lines 2+: SpecIdx1 SpecIdx2 NumShifts [sequence of space separated shift values]
%

fid = fopen(filename,'w');  if (fid<=0) fprintf(1,'ERROR: Unable to open %s!\n',filename); return; end;

numPairs = size(pairs,1);
fprintf(fid,'%d\n',numPairs);
for i=1:numPairs
    szList = size(shiftLists{i},1);
    fprintf(fid,'%d %d %d',double(pairs(i,1)),double(pairs(i,2)),szList);
    for j=1:szList
        fprintf(fid,' %.1f',shiftLists{i}(j));
    end
    fprintf(fid,'\n');
end

fclose(fid);
