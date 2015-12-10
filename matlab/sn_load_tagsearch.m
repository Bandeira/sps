function peptides = sn_load_tagsearch(specs, filename)
% function peptides = sn_load_tagsearch(specs, filename)

peptides = [];   lines = aux_load_lines(filename);   numLines = size(lines,1);   if numLines==0 return; end;
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

function lines = aux_load_lines(filename)
% function lines = aux_load_lines(filename)

fid = fopen(filename,'r');  if fid<=0 fprintf(1,'ERROR opening %s!',filename); lines={}; return; end;

buffer = fread(fid);
buffer = buffer(find(buffer~=13));  % Remove the annoying carriage return \r (DOS files)
lineBreaks = [0; find(buffer==10)];  if lineBreaks(size(lineBreaks,1))<size(buffer,1) lineBreaks=[lineBreaks; size(buffer,1)+1]; end;
lines = cell(size(lineBreaks,1)-1,1);
for lIdx = 1:size(lineBreaks,1)-1
    lines{lIdx} = sprintf('%s',char(buffer(lineBreaks(lIdx)+1:lineBreaks(lIdx+1)-1)));
end;

fclose(fid);
