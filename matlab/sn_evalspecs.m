function sn_evalspecs(contigs, specType, peakTol, mods, modMasses, addModMasses, numRegions, separateBY, curFNpref, excludeIso, maxIntensity)
% function sn_evalspecs(contigs, specType, peakTol, mods, modMasses, addModMasses, numRegions, separateBY, curFNpref, excludeIso, maxIntensity)

if isempty(curFNpref) 
    fprintf(1,'ERROR: Empty curFNpref! Why all the work if the results will not be saved?!'); return;
end;
if nargin==9 excludeIso=0; end;

if separateBY
	[bypercs,foo,blists,ipCountsB,ipCountsY] = sn_getbypercs(contigs, specType, peakTol, mods, modMasses, addModMasses, 0, 1, excludeIso, maxIntensity);
	idxB = find(sum(bypercs(:,5:6)')' >=  sum(bypercs(:,7:8)')');
	idxY = find(sum(bypercs(:,5:6)')' <  sum(bypercs(:,7:8)')');
    if numRegions>0
		[foo, regionProbsB]=sn_getbypercs(contigs(idxB,:), specType, peakTol, mods, modMasses, addModMasses, numRegions, 0, excludeIso, maxIntensity);
		[foo, regionProbsY]=sn_getbypercs(contigs(idxY,:), specType, peakTol, mods, modMasses, addModMasses, numRegions, 0, excludeIso, maxIntensity);
        eval(sprintf('%s_regionProbsB = regionProbsB; %s_regionProbsY = regionProbsY; %s_idxB = idxB; %s_idxY = idxY;',curFNpref,curFNpref,curFNpref,curFNpref));
%     else eval(sprintf('%s_idxB = idxB; %s_idxY = idxY;',curFNpref,curFNpref,curFNpref,curFNpref)); 
    end
else
    [bypercs, regionProbs,blists,ipCountsB,ipCountsY] = sn_getbypercs(contigs, specType, peakTol, mods, modMasses, addModMasses,numRegions, 1, excludeIso, maxIntensity);
%     eval(sprintf('%s_regionProbs = regionProbs;',curFNpref));
end;
[expInt,expPeaks]=aux_getsnr2(contigs, specType, mods, modMasses, addModMasses, peakTol);
% eval(sprintf('%s_bypercs = bypercs; %s_expInt = expInt; %s_expPeaks = expPeaks; %s_blists = blists; %s_ipCountsB = ipCountsB; %s_ipCountsY = ipCountsY;',curFNpref,curFNpref,curFNpref,curFNpref,curFNpref,curFNpref));

filename = sprintf('%s_sn_evalspecs',curFNpref);
if separateBY
    save(filename,'idxB','idxY','bypercs','expInt','expPeaks','blists','ipCountsB','ipCountsY');
else
    save(filename,'bypercs','expInt','expPeaks','blists','ipCountsB','ipCountsY');
end;

% Generate Excel .txt file
fid = fopen(sprintf('%s_sn_evalspecs_excel.txt',curFNpref),'w');
if fid<0 fprintf(1,'Error opening Excel file %s_sn_evalspecs_excel.txt',curFNpref);
else
    if separateBY
        fprintf(fid,'Explained intensity. Average main b/y = %.3f\n',mean([expInt(idxB,1);expInt(idxY,2)]));
        if isempty(idxB) t=0; else t=mean(expInt(idxB,:)); end; fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        if isempty(idxY) t=0; else t=mean(expInt(idxY,:)); end; fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        fprintf(fid,'Explained peaks. Average main b/y = %.3f\n',mean([expPeaks(idxB,1);expPeaks(idxY,2)]));
        if isempty(idxB) t=0; else t=mean(expPeaks(idxB,:)); end; fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        if isempty(idxY) t=0; else t=mean(expPeaks(idxY,:)); end; fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        fprintf(fid,'bypercs. Average main b/y = %.3f\n',mean([bypercs(idxB,1);bypercs(idxY,3)]));
        if isempty(idxB) t=0; else t=mean(bypercs(idxB,:)); end; fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        if isempty(idxY) t=0; else t=mean(bypercs(idxY,:)); end; fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        overallHist = sum([ipCountsB(idxB,:);ipCountsY(idxY,:)]);  if sum(overallHist)>0 overallHist = overallHist/sum(overallHist); end;
        fprintf(fid,'Inter-peak amino acid counts. Average main b/y = %.3f\n',sum(overallHist.*[1:6]));
        t=sum(ipCountsB(idxB,:)); if sum(t)>0 t=t/sum(t); end; fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        t=sum(ipCountsY(idxY,:)); if sum(t)>0 t=t/sum(t); end; fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        fprintf(fid,'regionProbsB\n');
        for i=1:numRegions fprintf(fid,'%.3f\t',regionProbsB(i,:));  fprintf(fid,'\n'); end;
        fprintf(fid,'regionProbsY\n');
        for i=1:numRegions fprintf(fid,'%.3f\t',regionProbsY(i,:));  fprintf(fid,'\n'); end;
    else
        fprintf(fid,'Explained intensity\n');
        t=mean(expInt); fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        fprintf(fid,'Explained peaks\n');
        t=mean(expPeaks); fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        fprintf(fid,'bypercs\n');
        t=mean(bypercs); fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        overallHist = sum([ipCountsB;ipCountsY]);   overallHist = overallHist/max(1,sum(overallHist));
        fprintf(fid,'Inter-peak amino acid counts. Average main b/y = %.3f\n',sum(overallHist.*[1:6]));
        t=sum(ipCountsB); t=t/max(1,sum(t)); fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        t=sum(ipCountsY); t=t/max(1,sum(t)); fprintf(fid,'%.3f\t',t); fprintf(fid,'\n');
        fprintf(fid,'regionProbs\n');
        for i=1:numRegions fprintf(fid,'%.3f\t',regionProbs(i,:)); fprintf(fid,'\n'); end;
    end
    
    fclose(fid);
end;

function [expInt, expPeaks, idxPeaks] = aux_getsnr2(specSet, specType, mods, modMasses, addModMasses, peakTol, idxPeaks, peptideMasses)
% function [expInt, expPeaks, idxPeaks] = aux_getsnr2(specSet, specType, mods, modMasses, addModMasses, peakTol, idxPeaks, peptideMasses)

global bProfile yProfile;

if nargin<8 peptideMasses=[]; end;
if size(specSet,2)==5 idxSpec=2; idxPM=3; idxPep=5; else idxSpec=3; idxPM=5; idxPep=7; end;
if strcmp(specType,'msms') peakOffset=0; else peakOffset=-1; end; 
ionsB = bProfile([8 11]',1)+peakOffset;    ionsBnl = bProfile([3:7]',1)+peakOffset;    ionsBz2 = bProfile([9 10]',1)+peakOffset;
ionsY = yProfile([6 8]',1)+peakOffset;     ionsYnl = yProfile([3:4]',1)+peakOffset;    ionsYz2 = yProfile([5 7]',1)+peakOffset;

numSpecs = size(specSet,1);   expInt = zeros(numSpecs,7);   expPeaks = zeros(numSpecs,8);   idxPeaks = cell(numSpecs,7);
for s=1:numSpecs
    if isempty(specSet{s,idxSpec}) | (isempty(specSet{s,idxPep}) & size(peptideMasses,1)<s) continue; end;
    if size(peptideMasses,1)>=s
        masses = unique([0 peptideMasses{s}(:,1)']);   
        masses = masses(2:length(masses))-masses(1:length(masses)-1);
    else 
        masses = sn_getmasses(specSet{s,idxPep},mods,modMasses,addModMasses); 
    end;
    szMasses = size(masses,2);
    idx = find(specSet{s,idxSpec}(:,1)>peakTol & specSet{s,idxSpec}(:,1)<sum(masses)-2*peakTol);    if isempty(idx) continue; end;
    specMasses = specSet{s,idxSpec}(idx,1)';   specScores = specSet{s,idxSpec}(idx,2)';    numPeaks = size(specMasses,2);   totIntensity = sum(specScores); if totIntensity==0 totIntensity=1; end;

    prefix = cumsum(masses(1:szMasses));    suffix = cumsum(masses(szMasses:-1:1));
    massesB = repmat(prefix,size(ionsB,1),1)+repmat(ionsB(:,1),1,szMasses);   massesB=reshape(massesB,size(ionsB,1)*szMasses,1);
    massesY = repmat(suffix,size(ionsY,1),1)+repmat(ionsY(:,1),1,szMasses);   massesY=reshape(massesY,size(ionsY,1)*szMasses,1);
    massesNL = [repmat(prefix,size(ionsBnl,1),1)+repmat(ionsBnl(:,1),1,szMasses); repmat(suffix,size(ionsYnl,1),1)+repmat(ionsYnl(:,1),1,szMasses)];
    massesNL = reshape(massesNL,size(massesNL,1)*szMasses,1);
    massesZ2 = [(repmat(prefix,size(ionsBz2,1),1)+repmat(ionsBz2(:,1),1,szMasses))/2; (repmat(suffix,size(ionsYz2,1),1)+repmat(ionsYz2(:,1),1,szMasses))/2];
    massesZ2 = reshape(massesZ2,size(massesZ2,1)*szMasses,1);
    
    % Masses of internal ions
    fragMasses = repmat(cumsum(masses),szMasses,1) - repmat(cumsum(masses)',1,szMasses);
    fragMasses = unique(fragMasses(find(fragMasses>0)))';   szFragMasses = size(fragMasses,2);
    if ~isempty(fragMasses)
        massesIfB = repmat(fragMasses,size(ionsB,1),1)+repmat(ionsB(:,1),1,szFragMasses);   massesIfB=reshape(massesIfB,size(ionsB,1)*szFragMasses,1);
        massesIfY = repmat(fragMasses,size(ionsY,1),1)+repmat(ionsY(:,1),1,szFragMasses);   massesIfY=reshape(massesIfY,size(ionsB,1)*szFragMasses,1);
        
        massesToTest = {massesB; massesY; massesNL; massesZ2; massesIfB; massesIfY};
    else massesToTest = {massesB; massesY; massesNL; massesZ2}; end;
    
    expPeaks(s,8) = numPeaks;   specPeaksIdx = [1:numPeaks];
    for m=1:size(massesToTest,1)
        if isempty(specMasses) break; end;
        diffs = min(abs(repmat(specMasses,size(massesToTest{m},1),1)-repmat(massesToTest{m},1,numPeaks)));
        match = find(diffs<=peakTol);   expInt(s,m)=sum(specScores(match));   expPeaks(s,m)=size(match,2);   idxPeaks{s,m}=specPeaksIdx(match);
        toKeep = find(diffs>peakTol);
        specMasses = specMasses(toKeep);   specScores = specScores(toKeep);    specPeaksIdx = specPeaksIdx(toKeep);
        numPeaks = size(toKeep,2);
    end
    expInt(s,7) = sum(specScores);
    expInt(s,:) = expInt(s,:) / totIntensity;
    expPeaks(s,7) = numPeaks;
    expPeaks(s,1:7) = expPeaks(s,1:7) / expPeaks(s,8);
    idxPeaks{s,7} = specPeaksIdx;
end
