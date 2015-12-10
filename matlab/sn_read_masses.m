function [res,masses,letters] = sn_read_masses(filename,masses,letters)
%  Reads amino acid masses from a text file of the format
%
%    X=<mass>
%
%  And outputs letters (char column array) and masses (column array of doubles).
%  Lines starting with '#' are ignored.
%
%  res = 1 if file loaded ok, 0 otherwise
%
fid = fopen(filename,'r'); if fid<=0 fprintf(1,'ERROR: Could not open %s!\n',filename); res=0; return; end; res=1;
masses = [];   letters = [];
while 1
    line = fgetl(fid);   szLine = length(line);
    if ~ischar(line) break; end;
    if szLine==0 | line(1)=='#' continue; end;
    v = str2num(line(3:szLine));
    if szLine<3 | line(2)~='=' | isempty(v) continue; end;
    
    letters = [letters; line(1)];    masses = [masses; v];
end

fclose(fid);
