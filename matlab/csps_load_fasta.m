function [seqs, status] = csps_load_fasta(filename, lineSep);
% function [seqs, status] = csps_load_fasta(filename, lineSep);
%
%  Returns the fasta sequences whose indices are in idx
%
%  Use lineSep = '\n' to keep the FASTA strings in the same format or '' to get one line sequences
%
%  seqs has 3 cols: protein identifier (col 1), sequence (col 2) and description (col 3)
%  status is set to -1 if there are errors processing the fasta file.
%

lines = sn_load_lines(filename);   numLines = size(lines,1);
if numLines==0 status=-1; seqs={}; return; end;

seqs = cell(numLines,3);

protIdx=0;
for lineIdx=1:numLines
    if mod(lineIdx,100000)==0 fprintf(1,'# lines = %d, # prots = %d\n',lineIdx, protIdx); end;
    line = lines{lineIdx};   if isempty(line) continue; end;
    if line(1)=='>'
        protIdx = protIdx + 1; 

        [dbId, line] = strtok(line(2:size(line,2)),'|');
        [protId, line] = strtok(line(2:size(line,2)),'|');
        seqs{protIdx,1} = sprintf('%s|%s',dbId, protId);
        seqs{protIdx,3} = line(2:size(line,2));
    else 
        if isempty(seqs{protIdx,2}) 
            seqs{protIdx,2}=line; 
        else
            seqs{protIdx,2} = sprintf('%s%s%s',seqs{protIdx,2},lineSep,line); 
        end;
    end
end;
seqs = seqs(1:protIdx,:);
status=0;
