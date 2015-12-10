function [peptidesDB, peptidesMatch, peptidesPos, results] = load_matchma(filename, separator, numEntries, numModStates)
% function [peptidesDB, peptidesMatch, peptidesPos, results] = load_matchma(filename, separator, numEntries, numModStates)
%
% Import matchma search results. Input file should be a matchma results file with _NO_ matched PRMs 
%
% peptides contains one line per entry and one col per number of mods with the highest scoring peptide in each case
%
% results is a cell variable with the full results: one entry per number of mods, each entry has as many
%   columns as there are columns in the matchma results file and contains the same data:
% Direction (b/y), Match score, Protein index ,Match start, Match end, Matched sequence, Matched peaks, Protein reference
%

NUM_HEADER_LINES = 2;
NUM_RESULT_COLS = 8;
PREFIX_SPEC = 'Spectrum ';
PREFIX_NUMMODS = 'Num mod/muts = ';
PREFIX_PRMS = 'PRMs interpreted as ';

fid=fopen(filename,'r'); if fid<=0 fprintf(1,'ERROR: Could not open %s!\n',filename); results={}; return; end;
for i=1:NUM_HEADER_LINES fgetl(fid); end;

results = cell(numEntries,numModStates);      peptidesPos = cell(1,numModStates);
peptidesDB = cell(numEntries,numModStates);   peptidesMatch = cell(numEntries,numModStates);
curEntry = 1;       str = fgetl(fid);
while(curEntry<=numEntries & (isempty(str) | str~=-1))
    while(~strncmp(str,PREFIX_NUMMODS,size(PREFIX_NUMMODS,2)) & (isempty(str) | str~=-1)) str = fgetl(fid); end;
    
    curMods=1;
    while(curMods<=numModStates & strncmp(str,PREFIX_NUMMODS,size(PREFIX_NUMMODS,2)))
        str = fgetl(fid);  % Get line after PREFIX_NUMMODS
        
        while(~strncmp(str,PREFIX_PRMS,size(PREFIX_PRMS,2)) & ~strncmp(str,PREFIX_NUMMODS,size(PREFIX_NUMMODS,2)) & ~strncmp(str,'Spectrum',8) & (isempty(str) | str~=-1)) str = fgetl(fid); end;
        firstResult = 1;
        while(strncmp(str,PREFIX_PRMS,size(PREFIX_PRMS,2)))
            curResult = cell(1,NUM_RESULT_COLS);
            [colStr, str] = strtok(str,separator);   
            curResult{1} = colStr(size(PREFIX_PRMS,2)+1:size(colStr,2));
            [colStr, str] = strtok(str,separator);   curResult{2} = str2num(colStr);
            [colStr, str] = strtok(str,separator);   curResult{3} = str2num(colStr);
            [colStr, str] = strtok(str,separator);   curResult{4} = str2num(colStr);
            [colStr, str] = strtok(str,separator);   curResult{5} = str2num(colStr);
            [colStr, str] = strtok(str,separator);   % Ignore length of matched sequence
            [colStr, str] = strtok(str,separator);   curResult{6} = colStr;
            [colStr, str] = strtok(str,separator);   curResult{7} = colStr;
            curResult{8} = str;
            results{curEntry,curMods}=[results{curEntry,curMods}; curResult];
            if firstResult 
                peptidesDB{curEntry,curMods}=curResult{6};    peptidesMatch{curEntry,curMods}=curResult{7}; 
                peptidesPos{curMods}=[peptidesPos{curMods}; [curResult{4} curResult{5}]];
                firstResult = 0;
            end;

            str = fgetl(fid);
        end
        curMods=curMods+1;
    end
    curEntry=curEntry+1;
end;

if curEntry<=numEntries fprintf(1,'ERROR: end of file reached before all expected entries were processed! Got %.0f out of %.0f expected (%s)!\n',curEntry,numEntries,filename); peptidesDB=[]; peptidesMatch=[]; peptidesPos=[]; results=[]; return; end; 

fclose(fid);