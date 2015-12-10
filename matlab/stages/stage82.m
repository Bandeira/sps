function stage82(contigs, specType, peakTol, mods, modMasses, addModMasses, numRegions, separateBY, curFNpref, excludeIso)
% function stage82(contigs, specType, peakTol, mods, modMasses, addModMasses, numRegions, separateBY, curFNpref, excludeIso)
%
%  Report diagnostic statistics on the starden alignment data
%
%  specType - 'prm' / 'msms'
%  mods, modMasses, addModMasses - for getmasses4
%  numRegions - for getBYpercs
%  separateBY - set to 1 if regionProbs are to be separated into b/y separate estimates, 0 otherwise
%  curFNpref  - filename/variables prefix
%  excludeIso - if excludeIso is specified and ==1 then isotopic peaks are not included in blists
%

if isempty(curFNpref) 
    fprintf(1,'ERROR: Empty curFNpref! Why all the work if the results will not be saved?!'); return;
end;
if nargin==9 excludeIso=0; end;

if separateBY
	[bypercs,foo,blists,ipCountsB,ipCountsY] = getBYPercs(contigs, specType, peakTol, mods, modMasses, addModMasses, 0, 1, excludeIso);
	idxB = find(sum(bypercs(:,5:6)') >= sum(bypercs(:,7:8)'));
	idxY = find(sum(bypercs(:,5:6)') <  sum(bypercs(:,7:8)'));
    if numRegions>0
		[foo, regionProbsB]=getBYPercs(contigs(idxB,:), specType, peakTol, mods, modMasses, addModMasses, numRegions, 0);
		[foo, regionProbsY]=getBYPercs(contigs(idxY,:), specType, peakTol, mods, modMasses, addModMasses, numRegions, 0);
        eval(sprintf('%s_regionProbsB = regionProbsB; %s_regionProbsY = regionProbsY; %s_idxB = idxB; %s_idxY = idxY;',curFNpref,curFNpref,curFNpref,curFNpref));
    else eval(sprintf('%s_idxB = idxB; %s_idxY = idxY;',curFNpref,curFNpref,curFNpref,curFNpref)); end
else
    [bypercs, regionProbs,blists,ipCountsB,ipCountsY] = getBYPercs(contigs, specType, peakTol, mods, modMasses, addModMasses,numRegions, 1, excludeIso);
    eval(sprintf('%s_regionProbs = regionProbs;',curFNpref));
end;
[expInt,expPeaks]=getsnr2(contigs, specType, mods, modMasses, addModMasses, peakTol);
eval(sprintf('%s_bypercs = bypercs; %s_expInt = expInt; %s_expPeaks = expPeaks; %s_blists = blists; %s_ipCountsB = ipCountsB; %s_ipCountsY = ipCountsY;',curFNpref,curFNpref,curFNpref,curFNpref,curFNpref,curFNpref));

v=version;   v = str2num(v(1));  if v>6 v=' -V6 '; else v=''; end;
eval(sprintf('save %s_stage82 %s %s_*;',curFNpref,v,curFNpref));

% Generate Excel .txt file
fid = fopen(sprintf('%s_stage82_excel.txt',curFNpref),'w');
if fid<0 fprintf(1,'Error opening Excel file %s_stage82_excel.txt',curFNpref);
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
