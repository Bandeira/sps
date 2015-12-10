function peptides = load_tagsearch(specs, filename)
%
%  Loads a set of peptides output by tagsearch
% 
%  peptides{i,1}{j,:} - Peptide (col.1), Tag score (col.2), Perc peptide score (col.3), ProtIdx (col.4), Protein (col.5)
%  peptides{i,2} - number of unique peptides
%  peptides{i,3} - unique peptides
%

peptides = [];   lines = load_lines(filename);   numLines = size(lines,1);   if numLines==0 return; end;
numSpecs = size(specs,1);   peptides = cell(numSpecs,3);
for lineIdx=2:numLines
    idxTab = find(lines{lineIdx}==char(9));  if length(idxTab)~=5 fprintf(1,'ERROR: Incorrect file format on line %d - too many TAB characters\n',lineIdx); peptides=[]; return; end;
    specIdx = str2num(lines{lineIdx}(1:idxTab(1)-1));    newPep = cell(1,5);    idxTab = [idxTab length(lines{lineIdx})];
    for idxField=1:5 newPep{idxField} = lines{lineIdx}(idxTab(idxField)+1:idxTab(idxField+1)-1); end
    for idxField=2:4 newPep{idxField} = str2num(newPep{idxField}); end;
    peptides{specIdx,1} = [peptides{specIdx,1}; newPep];
end

for specIdx=1:numSpecs
    if isempty(peptides{specIdx,1}) peptides{specIdx,2}=0; peptides{specIdx,3}=[]; 
    else peptides{specIdx,3} = unique(peptides{specIdx,1}(:,1)); peptides{specIdx,2}=size(peptides{specIdx,3},1); end;
end
