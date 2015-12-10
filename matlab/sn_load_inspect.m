function peptides = sn_load_inspect(prevPeps, filename, specFilenames, hasMSGF, loadAll)
% function peptides = sn_load_inspect(prevPeps, filename, specFilenames, hasMSGF, loadAll)
%
%  prevPeps - previous peptide annotations (only used to determine .dta filenames)
%  filename - filename of the inspect output file
%  specFilenames - filenames for spectra in specs, expected format is X:Y where X is a filename and Y is a 1-based spectrum index
%  hasMSGF - 0/1 depending on whether the results file contains the
%               additional MS-GF SpecProb field (default 0: no MS-GF probabilities)
%  loadAll - 0/1 indicator of whether to load all peptides per spectrum or
%               just the first one in the file (default 0: first peptide only)
%
%  peptides(i,:) - cell array with identified peptide (col.1), parent charge (col.2), indicator of whether an 
%                  annotation was provided for each spectrum (col.3), annotation MQScore (col.4), 
%                  annotation p-value (col.5), identified protein (col.6)
%
%  Expected columns:
%      SpectrumFile	- filename
%      Scan#	
%      Annotation	- peptide
%      Protein	
%      Charge	    - charge
%      MQScore	    - score
%      CutScore	
%      IntenseBY	
%      BYPresent	
%      Unused	
%      p-value	    - InsPecT p-value (FDR if after PValue.py)
%      DeltaCN	
%      DeltaCNOther	
%      RecordNumber	
%      DBFilePos	
%      SpecFilePos
%      SpecProb     - MG-GF p-value

if nargin<4 hasMSGF=0; loadAll=0; end;
if nargin<5 loadAll=0; end;

sep = char(9);  % TAB column separator
fid=fopen(filename,'r');  if fid<0 fprintf(1,'Error opening %s!\n',filename); return; end;

numSpecs=size(prevPeps,1);  peptides = cell(numSpecs,6+hasMSGF);    peptides(:,1) = prevPeps;
for i=1:numSpecs 
    peptides{i,6}='';   if isempty(peptides{i,1}) peptides{i,1}=''; end;
    if ~isempty(specFilenames)
        idx = find(specFilenames{i}=='/' | specFilenames{i}=='\');
        if ~isempty(idx) specFilenames{i} = specFilenames{i}(max(idx)+1:length(specFilenames{i})); end
    end;
end

line = fgetl(fid);  % Skip header line
line = fgetl(fid);   prevSpecname = '';   prevScanNum = -1;
while line~=-1 & ~isempty(line)
    if line(1)==sep  % Special processing for empty spectrum file name
        curSpecname = ''; str = line(2:length(line));
    else
        [tok,str] = strtok(line,sep);   curSpecname = tok(max([0 find(tok=='/' | tok=='\')])+1:size(tok,2));
    end;
    [tok,str] = strtok(str,sep); scanNum = str2num(tok);  % spectrum index in spectrum file
    if ~loadAll && strcmp(curSpecname,prevSpecname) && scanNum==prevScanNum line = fgetl(fid); continue; else prevSpecname=curSpecname; prevScanNum=scanNum; end;
    [tok,str] = strtok(str,sep); peptide = tok;  % peptide
    [tok,str] = strtok(str,sep); proteinID = tok;  % Protein ID
    [tok,str] = strtok(str,sep); charge = str2num(tok);  % Skip columns until charge

    if isempty(specFilenames) || isempty(specFilenames{i})
        idx = scanNum+1;
    else
        idx = min(find(strcmp(sprintf('%s:%d',curSpecname,scanNum+1),specFilenames)));  % +1 because InsPecT outputs 0-based spectrum indices
    end;
    if isempty(idx) 
        fprintf(1,'Warning: unable to find entry for %s (full line: %s)\n',curSpecname,line); line = fgetl(fid); continue; 
    end;
    idxPoints = find(peptide=='.');
    peptides{idx,1} = peptide(min(idxPoints)+1:max(idxPoints)-1);   peptides{idx,2} = charge;    peptides{idx,3} = 1;   peptides{idx,6} = proteinID;

	%  Convert InsPecT's CTG-57ALR+45TNV to CTG[-57]ALR[+45]TNV
    s = peptides{idx,1};  if isempty(s) continue; end;   
    s = strrep(s, 'phos', '+80');
    s=s(find(s<'a' | s>'z'));
    idx2 = find(s=='.'); if ~isempty(idx2) keep=[1 size(s,2)]; if min(idx2)<=2 keep(1)=min(idx2)+1; end; if max(idx2)>=size(s,2)-1 keep(2)=max(idx2)-1; end; s=s(keep(1):keep(2)); end;
    newS = ''; j=1;
    while j<=size(s,2)
        if s(j)>='A' & s(j)<='Z' newS=strcat(newS,s(j)); j=j+1;
        else
            newS = strcat(newS,'[');
            while j<=size(s,2) & (s(j)<'A' | s(j)>'Z') 
                if (s(j)=='-' | s(j)=='+') & newS(length(newS))~='[' newS=strcat(newS,']['); end;  % Two consecutive mods
                newS=strcat(newS,s(j)); j=j+1; 
            end;
            if j>size(s,2) newS=strcat(newS,']'); else newS = sprintf('%s]%c',newS,s(j)); j=j+1; end;
        end
    end
    peptides{idx,1}=newS;

    [tok,str] = strtok(str,sep);  peptides{idx,4} = str2num(tok);  % Get MQScore
    for i=1:5 [tok,str] = strtok(str,sep); end; peptides{idx,5} = str2num(tok);  % Get p-value
    
    if hasMSGF peptides{idx,7} = atof(str(max(find(str==sep))+1:length(str))); end
    
    line = fgetl(fid);
end
fclose(fid);
