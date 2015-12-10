function sn_save_lines(filename, lines, addLF)
% function sn_save_lines(filename, lines, addLF)
%  Writes a vector of lines to a file. Set addLF to 1 to add a \n at the end of every lines{i}.

fid = fopen(filename,'w');  if fid<=0 fprintf(1,'ERROR opening %s!',filename); return; end;

numLines = size(lines,1);   sizes = zeros(numLines,1);
for lIdx=1:numLines sizes(lIdx)=length(lines{lIdx})+addLF; end;

buffer = repmat(' ',1,sum(sizes));   sizes = [0; cumsum(sizes)];
for lIdx=1:numLines
    if addLF buffer(sizes(lIdx)+1:sizes(lIdx+1)-1) = lines{lIdx}; buffer(sizes(lIdx+1))=10; 
    else buffer(sizes(lIdx)+1:sizes(lIdx+1)) = lines{lIdx}; end
end;
fwrite(fid, buffer);

fclose(fid);
