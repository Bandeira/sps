function tags = load_tags(filename)
% function tags = load_tags(filename)
%
%  Loads tags generated with tags.cpp::ExtractTags()
%

lines = load_lines(filename);   ends = find(strcmp('END SPECTRUM',lines));   numSpecs = size(ends,1);
tags = cell(numSpecs,1);
for specIdx=1:numSpecs
    if specIdx==1 lineIdx=1; else lineIdx=ends(specIdx-1)+1; end;
    numTags = ends(specIdx)-1-lineIdx;   tags{specIdx}=cell(numTags,2);
    for tagIdx=1:numTags
        curStr = lines{lineIdx+tagIdx};
        split  = find(curStr==':');
        tags{specIdx}{tagIdx,1} = curStr(1:split-1);
        tags{specIdx}{tagIdx,2} = str2num(curStr(split+1:length(curStr)));
    end;
end;
