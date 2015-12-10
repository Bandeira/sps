function [expInt, expPeaks, idxPeaks] = getsnr2(specSet, specType, mods, modMasses, addModMasses, peakTol, idxPeaks, peptideMasses)
% function [expInt, expPeaks, idxPeaks] = getsnr2(specSet, specType, mods, modMasses, addModMasses, peakTol, idxPeaks, peptideMasses)
%
%  Like getsnr but returns explained intensity / # of peaks for
%    b/y/neutral losses/z=2/internal frags/other
%
%  mods, modMasses, addModMasses - passed directly to sn_getmasses
%  peptideMasses - if given then the masses in peptideMasses{i} are used instead of getmasses(peptide{i})
%  if idxPeaks==1 then it's converted to a cell structure containing the indices of the matched peaks per spectrum.
%       with as many cols as expInt (same meaning)
%
%  expInt, expPeaks - b/y/neutral losses/z=2/internal frags b/internal frags y/other : % exp. Intensity / % exp. peaks
%                     expPeaks has an additional column with # peaks in the spectrum
%
%  WARNING: Estimation of explained intensity / peaks is sequential (i.e. a peak is only tested as a y-ion if it was _not_ a b-ion)!
%

global bProfile yProfile;

if nargin<8 peptideMasses=[]; end;
if size(specSet,2)==5 idxSpec=2; idxPM=3; idxPep=5; else idxSpec=3; idxPM=5; idxPep=7; end;
% if strcmp(specType,'msms') peakOffset=0; else peakOffset=-1; end;
% ionsB = bProfile([8 11]',1)+peakOffset;    ionsBnl = bProfile([3:7]',1)+peakOffset;    ionsBz2 = bProfile([9 10]',1)+peakOffset;
% ionsY = yProfile([6 8]',1)+peakOffset;     ionsYnl = yProfile([2:3]',1)+peakOffset;    ionsYz2 = yProfile([5 7]',1)+peakOffset;
if strcmp(specType(1:2),'cz')
    if strcmp(specType,'czmsms') peakOffset=0; else peakOffset=-bProfile(12,1); end;
	ionsB = bProfile(12:13,1)+peakOffset;    ionsBnl = bProfile([3:7]',1)-bProfile(8,1)+bProfile(12,1)+peakOffset;    ionsBz2 = (bProfile([9 10]',1)-bProfile(8,1)+bProfile(12,1)+peakOffset)/2;
	ionsY = yProfile(9:10,1)+peakOffset;     ionsYnl = yProfile([3:4]',1)-yProfile(6,1)+yProfile(9,1)+peakOffset;     ionsYz2 = (yProfile([5 7]',1)-yProfile(6,1)+yProfile(9,1)+peakOffset)/2;
else
    if strcmp(specType,'msms') peakOffset=0; else peakOffset=-1; end;
	ionsB = bProfile([8 11]',1)+peakOffset;    ionsBnl = bProfile([3:7]',1)+peakOffset;    ionsBz2 = bProfile([9 10]',1)+peakOffset;
	ionsY = yProfile([6 8]',1)+peakOffset;     ionsYnl = yProfile([3:4]',1)+peakOffset;    ionsYz2 = yProfile([5 7]',1)+peakOffset;
end;

numSpecs = size(specSet,1);   expInt = zeros(numSpecs,7);   expPeaks = zeros(numSpecs,8);   idxPeaks = cell(numSpecs,7);
for s=1:numSpecs
    if isempty(specSet{s,idxSpec}) | (isempty(specSet{s,idxPep}) & size(peptideMasses,1)<s) continue; end;
    if size(peptideMasses,1)>=s
        masses = unique([0 peptideMasses{s}(:,1)']);
        masses = masses(2:length(masses))-masses(1:length(masses)-1);
%         prefix = peptideMasses{s};     suffix = unique(max(prefix)-[0 prefix]);     suffix=suffix(find(suffix>0));
%         szMasses = size(prefix,2);
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
%         if(m>1)  % Allow for peaks to be explained both as b AND y peaks
%             toKeep = find(diffs>peakTol);    if ~isempty(toRemove) toKeep=setdiff(toKeep,toRemove); toRemove=[]; end;
            toKeep = find(diffs>peakTol);
            specMasses = specMasses(toKeep);   specScores = specScores(toKeep);    specPeaksIdx = specPeaksIdx(toKeep);
            numPeaks = size(toKeep,2);
%         else toRemove = find(diffs<=peakTol); end;
    end
    expInt(s,7) = sum(specScores);
    expInt(s,:) = expInt(s,:) / totIntensity;
    expPeaks(s,7) = numPeaks;
    expPeaks(s,1:7) = expPeaks(s,1:7) / expPeaks(s,8);
    idxPeaks{s,7} = specPeaksIdx;
end