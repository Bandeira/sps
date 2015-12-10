function vSets = sn_report_specnets(specs, aligns, vSets, pp_matches, scores, threshExpInt, threshTP, graphvizCmd, labelpyDir, cysMass, filenames, proteins, specsMS, eiMSMS, reportFN)
% function vSets = sn_report_specnets(specs, aligns, vSets, pp_matches, scores, threshExpInt, threshTP, graphvizCmd, labelpyDir, cysMass, filenames, proteins, specsMS, eiMSMS, reportFN)

% specs     stars_dbannot       - annotated star spectra
% aligns    pairs_stars_annot   - pairs between annotated spectra (after propagation). Same as input to homglue?
%                                  used to define components, orient edges from smaller->larger parent masses and generate input files for graphviz
% pp_matches                    - pathproj matches: spectrum index (col.1) and propagation level (col.2)
% scores        eiPRM           - explained intensity / percent true b/y in PRM spectra
% threshExpInt  minEI           - minimum explained intensity threshold (to report an ID)
% threshTP      minTP           - minimum percent true b/y threshold (to report an ID)
% graphvizCmd   params.GRAPHVIZ_CMD - not used anymore
% labelpyDir    params.LABELPY_DIR - not used anymore
% cysMass                       - only needed for Label.py - not used anymore
% filenames                     - spectrum filenames to show on the report
% proteins                      - protein names for matched proteins
% specsMS   specs_raw           - only needed for Label.py - not used anymore
% eiMSMS    eiMSMS              - explained intensity / percent true b/y in MS/MS spectra 

if size(specs,2)==5 idxSpec=2; idxPM=3; idxPep=5; else idxSpec=3; idxPM=5; idxPep=7; end;
% fid=fopen(reportFN,'w');   if fid<=0 fprintf(1,'Error opening report_specnets.txt!\n'); return; end;

% Find all pairs of annotated spectra
idx_annot = [find(scores(:,1)>=threshExpInt & scores(:,2)>=threshTP)];
annotInfo = zeros(size(specs,1),3);    % Percent expInt (col.1), Percent TP (col.2), Annotation level (col.3)
annotInfo(:,1:2) = scores(:,1:2);    annotInfo(pp_matches(:,1),3) = pp_matches(:,2);
% idx_annot = [idx_annot; pp_matches(:,1)];  % Indices for initially/pathproj-annotated spectra

% Find all spectral networks (connected components)
pairs_annot = aligns(:,1:2); % pairs_annot = all spectral pairs
% [vSets,eSets] = aux_assembleE(double(pairs_annot));
numNets = size(vSets,1);   sizes = zeros(numNets,1); for i=1:numNets sizes(i)=size(vSets{i},2); end;

% Direct all edges from smaller to larger parent mass
for i=1:size(pairs_annot,1) if specs{pairs_annot(i,1),idxPM}>specs{pairs_annot(i,2),idxPM} pairs_annot(i,:) = pairs_annot(i,[2 1]); end; end;

% Output spectral networks report and graphs
dirContents = dir('.'); dirContents = {dirContents(:).name}'; 
if isempty(find(strcmp('graphs',dirContents))) status = dos('mkdir graphs'); end; 
if ~isempty(labelpyDir) labelpyUse=zeros(size(specsMS,1),1); end;
dirContents = dir('specnets'); dirContents = {dirContents(:).name}'; 
cd('graphs');   snIdx = zeros(size(specs,1),1);

report_lines = cell(sum(sizes)+1,1);   lineIdx=2;
report_lines{1} = sprintf('Network index\tSpectrum index\tSpectrum filename\tPeptide\tPropagation level\tPerc EI (score)\tPerc TP (score)\tPerc EI (MS/MS)\tProtein\n');
% fprintf(fid,'Network index\tSpectrum index\tSpectrum filename\tPeptide\tPropagation level\tPerc EI (score)\tPerc TP (score)\tPerc EI (MS/MS)\tProtein\n');
for i=1:numNets
    masses = [specs{vSets{i},idxPM}];   [foo,idxS]=sort(masses);   vSets{i}=vSets{i}(idxS);   
    vLabels = cell(sizes(i),3);   annotLevels = annotInfo(vSets{i},3)==0;
    if ~isempty(labelpyDir) & ~isempty(specsMS)
        if isempty(find(strcmp(sprintf('sn_%d',i),dirContents))) 
            status = mkdir(sprintf('../specnets/sn_%d',i));
        end;
    end;
    for j=1:sizes(i) 
        s = vSets{i}(j);   snIdx(s) = i;
        if ~isempty(filenames) filename = filenames{s}; else filename = ''; end;
        if ~isempty(proteins) protein = proteins{s}; else protein = ''; end;
        if ~isempty(eiMSMS) cur_eiMSMS = eiMSMS(s); else cur_eiMSMS = 0; end;
        if ~isempty(labelpyDir) labelpyUse(s)=1; end;
        report_lines{lineIdx} = sprintf('%d\t%d\t%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%s\n',i,s,filename,specs{s,5},annotInfo(s,3),100*annotInfo(s,1),100*annotInfo(s,2),100*cur_eiMSMS,protein);  lineIdx=lineIdx+1;
        vLabels{j,1} = s;   vLabels{j,2} = sprintf('%s',specs{s,5});   
        vLabels{j,3} = annotLevels(j);
    end;
    
%     aux_saveGraph(sprintf('specnet_%.0d.txt',i), pairs_annot(eSets{i},:), [], [], vLabels);
%     if ~isempty(graphvizCmd) status = dos(sprintf('%s -Tpng < specnet_%.0d.txt > specnet_%.0d.png',graphvizCmd,i,i)); if status<0 fprintf(1,'ERROR executing "%s"!\n',graphvizCmd); return; end; end;
end;
cd('..');
% fclose(fid);
sn_save_lines(reportFN,report_lines,0);


function [vSets, eSets] = aux_assembleE(aligns)
% function [vSets, eSets] = aux_assembleE(aligns)

numPeptides = size(unique(aligns(:,[1 2])),1);
vSets = cell(numPeptides,2);
contigsIdx = zeros(max(max(aligns(:,[1 2]))),1);
numContigSets = 0;

szAligns = size(aligns,1);
toProcess = ones(szAligns,1); idxAligns = [1:szAligns]';   tensLeft = floor(szAligns/10000);
curTime = 0; tic;
for idx=1:szAligns
    
    i = aligns(idx,1);
    j = aligns(idx,2);
    
    if contigsIdx(i)>0 & contigsIdx(j)>0 & contigsIdx(i)==contigsIdx(j) % ignore edges between vertices already in the same cluster - MAKES eSets INCOMPLETE !!!!!!!!
        continue; 
    end;  
    
    idxI = contigsIdx(i);
    idxJ = contigsIdx(j);
    
    if idxI==0 & idxJ==0 numContigSets=numContigSets+1; vSets{numContigSets,1} = [i j];   contigsIdx([i j])=numContigSets;   
    else if idxI==0 vSets{idxJ,1} = [vSets{idxJ,1} i];  contigsIdx(i)=idxJ;  
        else if idxJ==0 vSets{idxI,1} = [vSets{idxI,1} j];  contigsIdx(j)=idxI;
            else if idxI~=idxJ 
                    vSets{idxI,1} = unique([vSets{idxI,1} vSets{idxJ,1}]); 
                    contigsIdx(vSets{idxJ,1}) = idxI;
                    if idxJ<numContigSets
                        vSets{idxJ,1} = vSets{numContigSets,1}; 
                        contigsIdx(vSets{numContigSets,1}) = idxJ;
                    end
                    numContigSets=numContigSets-1;
                end;
            end;
        end;
    end;

    if max(contigsIdx)>numContigSets 
        a=1;
    end
    
    curTensLeft=floor(size(idxAligns,1)/10000); 
    if curTensLeft<tensLeft 
        t=toc; curTime=curTime+t; tic;
        fprintf(1,'%d aligns left, current 10000 took %.1f secs, ETA: %.1f secs...\n',size(idxAligns,1),t,curTensLeft*curTime/10000); tensLeft=curTensLeft; 
    end;
end; % while
vSets = vSets(1:numContigSets,:);
edgesSet = contigsIdx(aligns(:,1));   eSets = cell(size(vSets,1),1);
for i=1:size(vSets,1)  eSets{i} = find(edgesSet==i)'; end;

function idx = findContig(vertex, vSets)

idx=0;
for c=1:size(vSets,1)
    if ismember(vertex, vSets{c,1}) idx=c; return; end;
end;

function aux_saveGraph(filename, aligns, idxOk, idxSym, vLabels)
% function aux_saveGraph(filename, aligns, idxOk, idxSym, vLabels)

fid = fopen(filename,'w');
if fid<=0 fprintf(1,'Error opening %s\n',filename); return; end;

szAligns = size(aligns,1);
idxWr = setdiff([1:szAligns]', [idxOk; idxSym]);

fprintf(fid,'digraph G {');

labelsUsed = zeros(size(vLabels,1),1); labelsUsed(unique([aligns(:,1);aligns(:,2)]))=1;
for i=1:size(vLabels,1)  % used in sn_report_specnets for reporting spectral networks of modified peptides
    if ~labelsUsed(vLabels{i,1}) continue; end;
    if size(vLabels,2)==3 
        if vLabels{i,3}==1 fprintf(fid,'\t %d [ label ="%d: %s", color = black ]\n', vLabels{i,1}, vLabels{i,1}, vLabels{i,2});
        else fprintf(fid,'\t %d [ label ="%d: %s", color = red ]\n', vLabels{i,1}, vLabels{i,1}, vLabels{i,2}); end;
    else fprintf(fid,'\t %d [ label ="%d: %s" ]\n', vLabels{i,1}, vLabels{i,1}, vLabels{i,2}); end;
end

if size(aligns,2)>=3
	for i=1:size(idxOk,1) fprintf(fid,'\t %d -> %d [ color = green, style = bold, label = %.1f ];\n', aligns(idxOk(i),[1 2 3])); end;
	
	for i=1:size(idxSym,1) fprintf(fid,'\t %d -> %d [ color = red, style = bold, label = %.1f ];\n', aligns(idxSym(i),[1 2 3])); end;
	
	for i=1:size(idxWr,1) fprintf(fid,'\t %d -> %d [ color = blue, style = bold, label = %.1f ];\n', aligns(idxWr(i),[1 2 3])); end;
else
	for i=1:size(idxOk,1) fprintf(fid,'\t %d -> %d [ color = green, style = bold ];\n', aligns(idxOk(i),1:2)); end;
	for i=1:size(idxSym,1) fprintf(fid,'\t %d -> %d [ color = red, style = bold ];\n', aligns(idxSym(i),1:2)); end;
	for i=1:size(idxWr,1) fprintf(fid,'\t %d -> %d [ color = blue, style = bold ];\n', aligns(idxWr(i),1:2)); end;
end;

fprintf(fid,'}\n');
fclose(fid);

function peptides = aux_inspect_convertAnnotsRev(peptides,cysMass)
%  Converts annotations of the form AM[16]T[-18]R to Inspect's AM+16T-18R.

cysDiff = round(cysMass-160.0306482);
if cysDiff ~= 0
    if cysDiff<0 cysStr = sprintf('C%d',cysDiff); else cysStr = sprintf('C+%d',cysDiff); end;
else cysStr = ''; end;

numPeps = size(peptides,1);   
for p=1:numPeps
    if isempty(peptides{p}) continue; end;
    idxS=find(peptides{p}=='[');   idxE=find(peptides{p}==']');   numMods = size(idxS,2);   if numMods==0 continue; end;
    newPeptide = peptides{p}(1:idxS(1)-1);
    for m=1:numMods
        modMass = str2num(peptides{p}(idxS(m)+1:idxE(m)-1));
        if round(modMass)<0 newPeptide = sprintf('%s%d',newPeptide, round(modMass) );
        else newPeptide = sprintf('%s+%d',newPeptide, round(modMass) ); end;
        if m<numMods newPeptide = sprintf('%s%s',newPeptide, peptides{p}( idxE(m)+1:idxS(m+1)-1 ) ); end;  % Copy amino acids between annotations
    end
    peptides{p} = sprintf('%s%s',newPeptide, peptides{p}( idxE(numMods)+1 : size(peptides{p},2) ) );
    
    if peptides{p}(1)=='-' | peptides{p}(1)=='+'  % Move prefix +XX/-YY mods to after the first aa
        aaPos=min(find(peptides{p}>='A' & peptides{p}<'Z'));   aa = peptides{p}(aaPos);
        peptides{p}(2:aaPos) = peptides{p}(1:aaPos-1); peptides{p}(1) = aa;
    end
    
    if ~isempty(cysStr)  t = strrep(peptides{p},'C',cysStr);    peptides{p} = t;  end;
end;

function filenames = aux_savedta(data, scanNums, toDir)
%  Save a set of spectra to individual .dta files.

numSpecs = size(data,1);   filenames = cell(numSpecs,1);
for s=1:numSpecs
	if ~isempty(data{s,5}) filenames{s} = sprintf('%s/spec_%d_%s.dta', toDir, scanNums(s), data{s,5}); 
    else filenames{s} = sprintf('%s/spec_%d.dta', toDir, scanNums(s)); end; 

    fid = fopen(filenames{s},'w');    if fid<1 fprintf(1,'Error opening %s (index %d)!\n', filenames{s}, s); break; end;
    if ~isempty(data{s,4}) fprintf(fid,'%f %.0f\n', data{s,3}, data{s,4}); else fprintf(fid,'%f 2\n', data{s,3}); end;
    peaks = ''; 
    for p=1:size(data{s,2},1)
        t=sprintf('%s%f %f\n',peaks,data{s,2}(p,1:2));   peaks=t;
    end;
    fprintf(fid,peaks);
    fclose(fid);
end;
