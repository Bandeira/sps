function [countGIs, countProts] = fasta_countgis(filename, output)
% function [countGIs, countProts] = fasta_countgis(filename, output)

if nargin<2 output=1; end;

fid = fopen(filename,'r');   countGIs=-1;   countProts=-1;
if fid<0 fprintf(1,'ERROR: Cannot open file %s!\n', filename); return; end;

countGIs=0;  lines=0;   countProts=0;
line = fgetl(fid);  
while isempty(line) | line~=-1
    lines = lines+1; if output & mod(lines,10000)==0 fprintf(1,'# lines = %d, num. proteins = %d\n',lines, countProts); end;
    if size(line,2)>3 countGIs = countGIs + strcmp('gi', line(2:3)); countProts=countProts+strcmp(line(1),'>'); end;
    line = fgetl(fid);
end;
fclose(fid);
if output fprintf(1,'Total gi counts: %d\nTotal number of proteins: %d\n',countGIs,countProts); end;

