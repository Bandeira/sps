function [alignIdx, alignTxt, protsID] = clustalw_load_aln(filename)
% function [alignIdx, alignTxt, protsID] = clustalw_load_aln(filename)
%
%  Loads a .aln file as output by clustalw (v1.81) alignment of TWO proteins.
%
%  alignIdx    - 2-col vector of matched indices in prot 1/2 (cols 1/2), 0 if not matched
%  alignTxt{i} - alignment text for prots 1/2 (i=1/2) and match text (i=3)
%  protsID{i}  - protein IDs as in the .aln file (i=1/2)
%

alignIdx=[];   alignTxt={};   protsID={};
lines=sn_load_lines(filename);   numLines = size(lines,1);   if numLines<6 return; end;

alignTxt=cell(3,1);   protsID=cell(2,1);
[protsID{1}, alignTxt{1}] = strtok(lines{4},' ');    alignTxt{1} = alignTxt{1}(find(alignTxt{1}~=' '));
[protsID{2}, alignTxt{2}] = strtok(lines{5},' ');    alignTxt{2} = alignTxt{2}(find(alignTxt{2}~=' '));
mPos = findstr(lines{4},alignTxt{1});   alignTxt{3} = lines{6}(mPos:length(lines{6}));

% Read alignment text
lineIdx=8;
while lineIdx+2<=numLines
    [foo,txt] = strtok(lines{lineIdx},' ');      alignTxt{1} = strcat(alignTxt{1},txt(find(txt~=' ')));
    [foo,txt] = strtok(lines{lineIdx+1},' ');    alignTxt{2} = strcat(alignTxt{2},txt(find(txt~=' ')));
    alignTxt{3} = [alignTxt{3} lines{lineIdx+2}(mPos:length(lines{lineIdx+2}))];
    lineIdx = lineIdx+4;
end

% Get matched indices
alignIdx = zeros(length(alignTxt{1}),2);
idx = find(alignTxt{1}>='A' & alignTxt{1}<'Z');    alignIdx(idx,1) = [1:length(idx)]';
idx = find(alignTxt{2}>='A' & alignTxt{2}<'Z');    alignIdx(idx,2) = [1:length(idx)]';
