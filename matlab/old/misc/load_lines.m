function lines = load_lines(filename);
% function lines = load_lines(filename);
%
%  Reads everything from a text file and separates the lines into a cell column
%

fid = fopen(filename,'r');  if fid<=0 fprintf(1,'ERROR opening %s!',filename); lines={}; return; end;

buffer = fread(fid);
buffer = buffer(find(buffer~=13));  % Remove the annoying carriage return \r (DOS files)
lineBreaks = [0; find(buffer==10)];  if lineBreaks(size(lineBreaks,1))<size(buffer,1) lineBreaks=[lineBreaks; size(buffer,1)+1]; end;
lines = cell(size(lineBreaks,1)-1,1);
for lIdx = 1:size(lineBreaks,1)-1
    lines{lIdx} = sprintf('%s',char(buffer(lineBreaks(lIdx)+1:lineBreaks(lIdx+1)-1)));
end;

fclose(fid);
