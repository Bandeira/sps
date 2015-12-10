function peptides = inspect_loadAnnots3(prevPeps, filename, specFilenames, specFilenamesPrefix, fileType)
% function peptides = inspect_loadAnnots3(prevPeps, filename, specFilenames, specFilenamesPrefix, fileType)
%
%  Load inspect annotations in multiple-peptides-per-spectrum format. Similar to inspect_loadAnnots2 but does not require
%    the spectra as an input parameter.
%
%  prevPeps - previous peptide annotations (only used to determine .dta filenames)
%  filename - filename of the inspect output file
%  specFilenames       - one filename for each spectrum in specs (optional, set to [] for 'mgf' filetype)
%  specFilenamesPrefix - sometimes the spectrum filenames in inspect's output may be prepended with a common
%                          string. Common cases are when spectra were saved to dta files in a directory 'specs' or 'dtas' -
%                          in these cases set specFilenamesPrefix to 'specs/' or 'dtas/', respectively.
%  fileType - 'dta'/anything else influences usage of filenames ('dta') or filename+scan# (e.g. 'mgf') to match spectra to peptide IDs
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
%      p-value	    - p-value
%      DeltaCN	
%      DeltaCNOther	
%      RecordNumber	
%      DBFilePos	
%      SpecFilePos

%  spec_file_name, 8 other columns, spec_annotation, 2 other columns
%
%
%  --- >>> For the simpler 5 column format: 
%              spec_file_name, peptideID, score, p-value, peptide charge
%  --- >>> Just use data_annotations=read_csv(filename,separator,numLines,5);  <<< ---
%

if nargin<5 fileType='dta'; end;

sep = char(9);  % TAB column separator

fid=fopen(filename,'r');  if fid<0 fprintf(1,'Error opening %s!\n',filename); return; end;

numSpecs=size(prevPeps,1);  peptides = cell(numSpecs,6);    peptides(:,1) = prevPeps;   for i=1:numSpecs if isempty(peptides{i,1}) peptides{i,1}=''; end; peptides{i,6}=''; end;
if strcmp(fileType,'dta') & isempty(specFilenames)
    specFilenames = cell(numSpecs,1);
    for i=1:numSpecs if ~isempty(prevPeps{i}) specFilenames{i} = sprintf('spec_%d_%s.dta',i,prevPeps{i}); else specFilenames{i} = sprintf('spec_%d.dta',i); end; end;
end;
if ~isempty(specFilenamesPrefix)
    for i=1:numSpecs specFilenames{i} = sprintf('%s%s',specFilenamesPrefix,specFilenames{i}); end; 
end;

line = fgetl(fid);   prevSpecname = '';   prevScanNum = -1;
while line~=-1 & ~isempty(line)
    [tok,str] = strtok(line,sep);   curSpecname = tok(max([1 find(tok=='/' | tok=='\')])+1:size(tok,2));
%     for i=1:9 [tok,str] = strtok(str,sep); end;  peptide =tok;  % Skip columns until annotated peptide
    [tok,str] = strtok(str,sep); scanNum = str2num(tok);  % spectrum index (inside spectrum file, always ==1 for dta files)
    if strcmp(curSpecname,prevSpecname) & scanNum==prevScanNum line = fgetl(fid); continue; else prevSpecname=curSpecname; prevScanNum=scanNum; end;
    [tok,str] = strtok(str,sep); peptide = tok;  % peptide
    [tok,str] = strtok(str,sep); proteinID = tok;  % Protein ID
    [tok,str] = strtok(str,sep); charge = str2num(tok);  % Skip columns until charge

    if strcmp(fileType,'dta') idx = min(find(strcmp(curSpecname,specFilenames))); else idx=scanNum+1; end;
    if isempty(idx) 
        fprintf(1,'Warning: unable to find entry for %s (full line: %s)\n',curSpecname,line); line = fgetl(fid); continue; 
    end;
    idxPoints = find(peptide=='.');
    peptides{idx,1} = peptide(min(idxPoints)+1:max(idxPoints)-1);   peptides{idx,2} = charge;    peptides{idx,3} = 1;   peptides{idx,6} = proteinID;

    [tok,str] = strtok(str,sep);  peptides{idx,4} = str2num(tok);  % Get MQScore
    for i=1:5 [tok,str] = strtok(str,sep); end; peptides{idx,5} = str2num(tok);  % Get p-value
    
    line = fgetl(fid);
end
fclose(fid);
fprintf(1,'Finishing with %d spectra left unchanged\n',numSpecs-length(find([peptides{:,3}]==1)));
