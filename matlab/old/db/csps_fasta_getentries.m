function [seqs, status] = fasta_getentries(idx, filename, lineSep);
% function [seqs, status] = fasta_getentries(idx, filename, lineSep);
%
%  Returns the fasta sequences whose indices are in idx
%
%  Use lineSep = '\n' to keep the FASTA strings in the same format or '' to get one line sequences
%
%  seqs has 3 cols: protein identifier (col 1), sequence (col 2) and description (col 3)
%  status is set to -1 if there are errors processing the fasta file.
%
idx = sort(idx);

fid = fopen(filename,'r'); seqs={}; status=-1;
if fid<0 fprintf(1,'Cannot open file %s!\n', filename); return; end;

seqs = cell(size(idx,1),3);

count=0;  lines=0;   idxI=1;
line = fgetl(fid);  
while line~=-1
    lines = lines+1; if mod(lines,100000)==0 fprintf(1,'# lines = %d, count = %d, idx left = %d\n',lines, count, size(find(idx>count),1)); end;
    if line(1)=='>'
        count = count + 1; 
        if count>max(idx) return; end;
        if count>idx(idxI) idxI = idxI+1; end;

        [dbId, line] = strtok(line(2:size(line,2)),'|');
        [protId, line] = strtok(line(2:size(line,2)),'|');
        if count==idx(idxI)
            seqs{idxI,1} = sprintf('%s|%s',dbId, protId);
            seqs{idxI,3} = line(2:size(line,2));
        end
    else 
        if count==idx(idxI) if isempty(seqs{idxI,2}) seqs{idxI,2}=line; else seqs{idxI,2} = sprintf('%s%s%s',seqs{idxI,2},lineSep,line); end; end;
    end
    line = fgetl(fid); while isempty(line) line = fgetl(fid); end;
end;

fclose(fid);
status=0;
