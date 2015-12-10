function contents = read_csv(filename,separator, numCols)
% function contents = read_csv(filename,separator, numCols)
%
%  Read the contents of a csv file into a 2D cell structure with string entries
%

lines=load_lines(filename);   numLines=size(lines,1);   if numLines==0 contents={}; return; end;

if strcmp(separator,'\t')==1 
    separator=char(9); 
end;
if nargin<3 
    numCols=length(find(lines{1}==separator))+1; 
end;

contents = cell(numLines,numCols);
for i=1:numLines
    if ~isempty(lines{i})
        idx = [0 find(lines{i}==separator) size(lines{i},2)+1];
        curLim = min(numCols,size(idx,2)-1);
        for j=1:curLim
            contents{i,j}=lines{i}(idx(j)+1:idx(j+1)-1);
        end
        for j=curLim+1:numCols contents{i,j} = ''; end;  % Fill missing values with empty strings
    end
end;
