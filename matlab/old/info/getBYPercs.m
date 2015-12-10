function [percs, regionProbs, blists, ipCountsB, ipCountsY] = getBYPercs(specSet, specType, peakTol, mods, modMasses, addModMasses, numRegions, getBLists, excludeIso, maxIntensity)
% function [percs, regionProbs, blists, ipCountsB, ipCountsY] = getBYPercs(specSet, specType, peakTol, mods, modMasses, addModMasses, numRegions, getBLists, excludeIso, maxIntensity)
%
%  Find what percentages of true b/y ions are represented in each spectrum
%
%  percs(i,:) - [% true b, % true b-iso, % true y, % true y-iso, % explained intensity for b/bIso/y/yIso (4 cols), # of b/y peaks considered (col.9)]
%  blists{i}  - Contains a list of (b-ion index, % intensitity) pairs containing all present b-ions
%  regionProbs - matrix numRegions lines by 4 cols: [p(b) p(bIso) p(y) p(yIso)] per region
%  ipCountsB(i,j) - Number of times that consecutive matched b peaks differ by j amino acids (j<=6)
%  ipCountsY(i,j) - Number of times that consecutive matched y peaks differ by j amino acids (j<=6)
%  excludeIso - if excludeIso is specified and ==1 then isotopic peaks are not included in blists
%  maxIntensity   - if maxIntensity is specified and ==1 then each theoretical peak's %explained intensity is the maximum of all peaks within tolerance (instead of the usual sum)
%

global bProfile yProfile;
if nargin<9  excludeIso=0; end;
if nargin<10 maxIntensity=0; end;

numSpecs = size(specSet,1);   percs = zeros(numSpecs,2);
if size(specSet,2)==5 specIdx=2; pepIdx=5; else specIdx=3; pepIdx=7; end;
% if strcmp(specType,'msms') peaksOffset=1; else peaksOffset=0; end;
% if strcmp(specType,'msms') peaksOffset=0; else peaksOffset=1.0072763; end;
if strcmp(specType,'czprm') peaksOffsetB=0; peaksOffsetY=yProfile(9,1)-bProfile(12,1); else 
    if strcmp(specType,'czmsms') peaksOffsetB=bProfile(12,1); peaksOffsetY=yProfile(9,1); else
        if strcmp(specType,'msms') peaksOffsetB=bProfile(8,1); peaksOffsetY=yProfile(6,1); else peaksOffsetB=0; peaksOffsetY=yProfile(6,1)-1.0072763; end;
    end
end;

percs = zeros(numSpecs,9);   regionProbs = zeros(numRegions,4);   if ~isempty(getBLists) & getBLists==1 blists=cell(numSpecs,1); else getBLists=0; end;
ipCountsB = zeros(numSpecs,6);   ipCountsY = zeros(numSpecs,6);

for s=1:numSpecs
    if isempty(specSet{s,pepIdx}) | isempty(specSet{s,specIdx}) continue; end;
    totalIntensity = sum(specSet{s,specIdx}(find(specSet{s,specIdx}(:,1)>5),2));   if totalIntensity==0 totalIntensity=1; end;
    masses = getmasses4(specSet{s,pepIdx},mods,modMasses,addModMasses);   szMasses = size(masses,2);   percs(s,9) = szMasses-1;   if szMasses==1 continue; end; % bn/yn are always trivially matched
    specSet{s,specIdx} = double(specSet{s,specIdx}(find(specSet{s,specIdx}(:,1)>20 & specSet{s,specIdx}(:,1)<=sum(masses)+peakTol),:));  % Ignore b0/y0/bn/yn and peaks with masses larger than the peptide mass
    massesB = cumsum(masses(1:szMasses-1))+peaksOffsetB;   massesBiso = massesB+bProfile(11,1)-bProfile(8,1);           % Also look for isotopic peaks
    massesY = cumsum(masses(szMasses:-1:2))+peaksOffsetY;  massesYiso = massesY+yProfile(8,1)-yProfile(6,1);  
    
    repSpec = repmat(specSet{s,specIdx}(:,1),1,szMasses-1);   if getBLists==1 | maxIntensity==1 repScores = repmat(specSet{s,specIdx}(:,2),1,szMasses-1); end;
    
    szSpec = size(specSet{s,specIdx},1);   idxSpec = reshape([1:szSpec;1:szSpec],2*szSpec,1);
    diffs = abs(repmat(massesB,szSpec,1) - repSpec);
    if size(diffs,1)>1 matchesB = min(diffs) <= peakTol; else matchesB = diffs <= peakTol; end;   matchesSpecIdx = find(min(diffs')<=peakTol);
    matchedIdx=find(matchesB==1);   ipCountsB(s,:)=hist(matchedIdx(2:size(matchedIdx,2))-matchedIdx(1:size(matchedIdx,2)-1),[1:6]);
    if(maxIntensity) percs(s,[1 5]) = [size(find(matchesB==1),2)/(szMasses-1) sum(max((diffs<=peakTol).*repScores))];
    else percs(s,[1 5]) = [size(find(matchesB==1),2)/(szMasses-1) sum(specSet{s,specIdx}(matchesSpecIdx,2))]; end;
    if getBLists==1 & szSpec>1 blists{s} = max((diffs<=peakTol).*repScores); end;
    
    diffs = abs(repmat(massesBiso,szSpec,1) - repSpec);
    if size(diffs,1)>1 matchesBiso = min(diffs) <= peakTol; else matchesBiso = diffs <= peakTol; end;  newMatches = find(min(diffs')<=peakTol);
%    percs(s,[2 6]) = [size(find(matchesBiso==1),2)/(szMasses-1) sum(specSet{s,specIdx}(setdiff(newMatches,matchesSpecIdx),2))];
    if(maxIntensity) percs(s,[2 6]) = [size(find(matchesBiso==1),2)/(szMasses-1) sum(max((diffs<=peakTol).*repScores))];
    else percs(s,[2 6]) = [size(find(matchesBiso==1),2)/(szMasses-1) sum(specSet{s,specIdx}(newMatches,2))]; end;
    matchesSpecIdx = [matchesSpecIdx newMatches];
    if getBLists==1 & szSpec>1 
        if (excludeIso==1)  blists{s} = [matchesB'.*[1:szMasses-1]' blists{s}']; 
        else blists{s} = [(matchesBiso|matchesB)'.*[1:szMasses-1]' max([(diffs<=peakTol).*repScores;blists{s}])']; end;
        blists{s}=blists{s}(find(blists{s}(:,1)>0),:); blists{s}(:,2)=blists{s}(:,2)/totalIntensity;
    end;
    
    diffs = abs(repmat(massesY,szSpec,1) - repSpec);
    if size(diffs,1)>1 matchesY = min(diffs) <= peakTol; else matchesY = diffs <= peakTol; end;  newMatches = find(min(diffs')<=peakTol);
    matchedIdx=find(matchesY==1);   ipCountsY(s,:)=hist(matchedIdx(2:size(matchedIdx,2))-matchedIdx(1:size(matchedIdx,2)-1),[1:6]);
%    percs(s,[3 7]) = [size(find(matchesY==1),2)/(szMasses-1) sum(specSet{s,specIdx}(setdiff(newMatches,matchesSpecIdx),2))];
    if(maxIntensity) percs(s,[3 7]) = [size(find(matchesY==1),2)/(szMasses-1) sum(max((diffs<=peakTol).*repScores))];
    else percs(s,[3 7]) = [size(find(matchesY==1),2)/(szMasses-1) sum(specSet{s,specIdx}(newMatches,2))]; end;
    matchesSpecIdx = [matchesSpecIdx newMatches];
    
    diffs = abs(repmat(massesYiso,szSpec,1) - repSpec);
    if size(diffs,1)>1 matchesYiso = min(diffs) <= peakTol; else matchesYiso = diffs <= peakTol; end;  newMatches = find(min(diffs')<=peakTol);
%    percs(s,[4 8]) = [size(find(matchesYiso==1),2)/(szMasses-1) sum(specSet{s,specIdx}(setdiff(newMatches,matchesSpecIdx),2))];
    if(maxIntensity) percs(s,[4 8]) = [size(find(matchesYiso==1),2)/(szMasses-1) sum(max((diffs<=peakTol).*repScores))];
    else percs(s,[4 8]) = [size(find(matchesYiso==1),2)/(szMasses-1) sum(specSet{s,specIdx}(newMatches,2))]; end;
    matchesSpecIdx = [matchesSpecIdx newMatches];
    
%     diffs = min(abs(repmat(massesY,szSpec,1) - repmat(specSet{s,specIdx}(idxSpec,1),1,szMasses-1))');
%     percs(s,3:4) = size(find(diffs<=peakTol),2)/(szMasses-1);
    
    percs(s,5:8) = percs(s,5:8)./totalIntensity;
%    percs(s,5:8) = percs(s,5:8)./sum(specSet{s,specIdx}(:,2));
    
    if numRegions==0 continue; end;
    regionSize = (sum(masses)+2*peakTol)/numRegions;   regionCounts = zeros(numRegions,4);
    massesAll = [ceil(massesB./regionSize)' ceil(massesBiso./regionSize)' ceil(massesY./regionSize)' ceil(massesYiso./regionSize)'];
    massesAll(find(massesAll>numRegions)) = numRegions;   massesAll(find(massesAll<1)) = 1;
    matchesAll= [matchesB' matchesBiso' matchesY' matchesYiso'];
    for i=1:4
        for j=1:szMasses-1
            regionCounts(massesAll(j,i),i) = regionCounts(massesAll(j,i),i)+matchesAll(j,i);
        end
    end
    regionMax = hist(massesAll,[1:numRegions]);   idxZeros = find(regionMax==0);   regionMax(idxZeros)=1;
    regionCounts = regionCounts./regionMax;       regionCounts(idxZeros)=1;
    regionProbs = regionProbs + regionCounts;
end;
if numSpecs>0 regionProbs = regionProbs./numSpecs; end;
