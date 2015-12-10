function sn_genstars(contigs, peakTol, pmTol, pValue, noiseProb, scoreMergeType, alignedPairs, pairspectraFN, modsFN, specsFNprefix)
% function sn_genstars(contigs, peakTol, pmTol, pValue, noiseProb, scoreMergeType, alignedPairs, pairspectraFN, modsFN, specsFNprefix)

% Load star spectra and mod positions
modPositions = sn_load_binarray(modsFN, 'float');
starden_contigs = sn_load_pklbin(pairspectraFN);

contigs_idx = reshape(alignedPairs',2*size(alignedPairs,1),1)';   contigs_idx_other = reshape(alignedPairs(:,[2 1])',2*size(alignedPairs,1),1)';
unique_contigs = unique(contigs_idx);   unique_contigs_size = size(unique_contigs,2);   vSets = cell(unique_contigs_size,1);
sn_save_binarray('stars_indices.bin', unique_contigs, 'int32');
for i=1:unique_contigs_size
    vSets{i} = find(contigs_idx==unique_contigs(i));
end

% Copy peptide annotations (if existing)
if size(contigs,2)==5 idxPM=3; idxPep=5; else idxPM=5; idxPep=7; end;
starden_contigs(:,5) = contigs(contigs_idx,idxPep);

masses = [contigs{:,idxPM}]';  masses = masses(alignedPairs');
otherPM = reshape([masses(2,:); masses(1,:)],size(contigs_idx,2),1)';
[dirs, consensusPA, endpoints] = aux_orientPairedSpectra(starden_contigs, vSets, modPositions, otherPM, peakTol, pmTol, contigs(contigs_idx,idxPep), pValue, noiseProb, 2, scoreMergeType);
sn_save_pklbin(consensusPA,sprintf('%sstars_only.pklbin',specsFNprefix));

% Create consensusPAext to include consensusPA and non-deconvolved spectra while keeping every spectrum's original index
if size(contigs,2)==5 consensusPAext = contigs; else consensusPAext = contigs(:,[1 3 5:7]); end;
consensusPAext(unique_contigs,:) = consensusPA;
sn_save_pklbin(consensusPAext,sprintf('%sstars.pklbin',specsFNprefix));

function [dirs, consensus, endpoints] = aux_orientPairedSpectra(specSet, vSets, modOffsets, otherPM, peakTol, pmTol, peptides, pValue, probNoise, mergeType, scoreMergeType, resolution)
% function [dirs, consensus, endpoints] = aux_orientPairedSpectra(specSet, vSets, modOffsets, otherPM, peakTol, pmTol, peptides, pValue, probNoise, mergeType, scoreMergeType, resolution)

if nargin<12 resolution=0.1; end;

numSets = size(vSets,1);   dirs = cell(numSets,1);   consensus = cell(numSets,5);   endpoints = cell(numSets,1);
if resolution==0.1 pmTolRange = [-pmTol:.5:pmTol]; else pmTolRange = [-pmTol:resolution:pmTol]; end;
szRange = size(pmTolRange,2);
for i=1:numSets
    if isempty(vSets{i}) continue; end;
    specs = specSet(vSets{i},:);   numSpecs = size(specs,1);   pmMass = specs{1,3}-1;    dirs{i}=zeros(numSpecs,1);
    endpoints{i} = [];  % Separate structure to hold endpoint PRMs until de-novo interpretation. Can't tell between b/y internal/endpoint PRMs at this point
                        %  All endpoints are stored assuming prefix PRMs are getting matched between the 2 spectra
    if ~isempty(modOffsets) curOffsets = modOffsets(ceil(vSets{i}/2),:); end;  % division by 2 converts from specIdx to pairIdx

    baseSpec = [specs(1,:);specs(1,:)];  % Second entry will contain the reversed spectrum (line below)
    baseSpec{2,2}(:,1) = pmMass-baseSpec{2,2}(:,1);   [foo, idxS]=sort(baseSpec{2,2}(:,1));   baseSpec{2,2}=baseSpec{2,2}(idxS,:);

    for s=1:numSpecs
        if isempty(specs{s,2}) continue; end;
        curEndpoints = [];
        if ~isempty(modOffsets) & curOffsets(s)>=0 & otherPM(vSets{i}(s))<pmMass+1-pmTol  % second condition because only the larger spectrum can gain endpoint peaks
            if curOffsets(s)<=pmTol curEndpoints = [pmMass+1-otherPM(vSets{i}(s)) max(specs{s,2}(:,2))]; end;    % Mod at the start
            if curOffsets(s)>=otherPM(vSets{i}(s))-pmTol curEndpoints = [otherPM(vSets{i}(s))-19 max(specs{s,2}(:,2))]; end;  % Mod at the end
        end;
        if s==1 endpoints{i} = curEndpoints; continue; end % this line added to allow the endpoints code above to be run for s=1 also

        [idx1, idx2] = aux_findMatchPeaks(baseSpec{1,2}, specs{s,2}, 0, peakTol);
        scoreB = sum(baseSpec{1,2}(idx1,2)) + sum(specs{s,2}(idx2,2));

        scoreY=0;  bestShift = pmTol;
        for r=1:szRange
            [idx1, idx2] = aux_findMatchPeaks(baseSpec{2,2}, specs{s,2}, pmTolRange(r), peakTol);
            curScoreY = sum(baseSpec{2,2}(idx1,2)) + sum(specs{s,2}(idx2,2));
            if curScoreY==scoreY & abs(pmTolRange(r))<bestShift bestShift=pmTolRange(r); end
            if curScoreY>scoreY scoreY=curScoreY; bestShift=pmTolRange(r); end;
        end

        if (scoreY>scoreB)
            specs{s,2}(:,1) = pmMass-specs{s,2}(:,1);   [foo, idxS]=sort(specs{s,2}(:,1));   specs{s,2}=specs{s,2}(idxS,:);
            if ~isempty(curEndpoints) curEndpoints(1,1) = pmMass-18-curEndpoints(1,1); end;  % At most one endpoint per spectrum
        end;
        endpoints{i} = [endpoints{i}; curEndpoints];
        dirs{i}(s) = scoreY>scoreB;
    end;

    consensus(i,:) = aux_getClusterConsensus([specs(:,1:4) peptides(vSets{i})], {1:size(specs,1)}, peakTol, pValue, probNoise, mergeType, scoreMergeType, resolution);
end;

function [idx1, idx2] = aux_findMatchPeaks(prmList1, prmList2, shift, tolerance)
% function [idx1, idx2] = aux_findMatchPeaks(prmList1, prmList2, shift, tolerance)

prmList1 = [0 0; prmList1];         prmList2 = [0 0; prmList2(:,1)+shift prmList2(:,2)];
szList1 = size(prmList1,1);         szList2 = size(prmList2,1);
scores = zeros(szList1, szList2);   preds = cell(szList1, szList2);   preds{1,1} = [1 1; 0 0]; % Best tuple so far, predecessor

pairs = zeros(szList1*szList2,9);  % oversized list of pairs: prm i, prm j, bestScoreUntilIJ, bestPairUntilIJ (2 cols), predecessors (2 cols), predIdx, matchLength
pairs(1,:) = [1 1 0 0 0 0 0 0 0];
curPair = 2;

startJ = 2;  % Prm at which to start looking for PRMs (in prmList2) matching prmList1(i,1)
for i=2:szList1
    for j=startJ:szList2
        if (prmList1(i,1)>prmList2(j,1)+tolerance) startJ = j+1; continue; end;
        if (prmList1(i,1)<prmList2(j,1)-tolerance) break; end;

        prevPair=curPair-1; while pairs(prevPair,1)>=i | pairs(prevPair,2)>=j prevPair=prevPair-1; end;  % find (i-1,j-1)

        if pairs(prevPair,3)+prmList1(i,2)+prmList2(j,2)>pairs(curPair-1,3)
            pairs(curPair,:) = [i j pairs(prevPair,3)+prmList1(i,2)+prmList2(j,2) i j pairs(prevPair,[4 5]) prevPair pairs(prevPair,9)+1];
        else
            pairs(curPair,:) = [i j pairs(curPair-1,3:9)];
        end;

        curPair = curPair + 1;
    end
end

curPair = curPair-1;      szKeep = pairs(curPair,9);      toKeep = zeros(szKeep,2);
for i=szKeep:-1:1
    toKeep(i,:) = pairs(curPair,[4 5]);   curPair = pairs(curPair,8);
end

idx1 = toKeep(:,1)-1;
idx2 = toKeep(:,2)-1;

function consensus = aux_getClusterConsensus(specSet, vSets, windowSize, pValue, noiseProb, mergeType, scoreMergeType, resolution)
% function consensus = aux_getClusterConsensus(specSet, vSets, windowSize, pValue, noiseProb, mergeType, scoreMergeType, resolution)

if nargin<8 resolution=0.1; end;

if size(specSet,2)==5 idxSpec=2; idxPM=3; idxPep=5; isContig=0; numCols=5; else idxSpec=3; idxPM=5; idxPep=7; isContig=1; numCols=7; end;
numClusters = size(vSets,1);    consensus = cell(numClusters,numCols);
numBinsInWindow = 1+2*windowSize/resolution;
numRegions = size(noiseProb,2); pSuccess = 1-(1-noiseProb).^numBinsInWindow;   pValue = pValue * ones(1,numRegions);
for i=1:numClusters
    regions = round(cumsum(ones(1,numRegions)*max([specSet{vSets{i,1},idxPM}])/numRegions)/resolution);

    toKeep = binoinv( pValue, size(vSets{i,1},2)*ones(1,numRegions) , pSuccess);  % toKeep(i)-minimum number of matches in the i-th region
    [counts,matches]=aux_ciPeakCoherence(specSet(vSets{i,1},:),isContig,windowSize,0,resolution);

    consensus{i,idxSpec} = [];   regionStart = 0;
    for j=1:numRegions
        idx = regionStart+find(matches((regionStart+1):min(regions(j),size(matches,1)),1)>=max(1,toKeep(j)));
        if strcmp(scoreMergeType,'sum')
            consensus{i,idxSpec} = [consensus{i,idxSpec}; [idx*resolution matches(idx,3)]];
        else
            consensus{i,idxSpec} = [consensus{i,idxSpec}; [idx*resolution matches(idx,3)./matches(idx,1)]];
        end
        regionStart = regions(j);
    end
    if mergeType==1 consensus{i,idxSpec} = aux_weightedMerge(consensus{i,idxSpec}, windowSize, 'max'); end;
    if mergeType==2
        if ~isempty(consensus{i,idxSpec}) consensus{i,idxSpec} = aux_mergeConsensusPeaks(consensus{i,idxSpec},round(windowSize/resolution),resolution); end;
    end;

    consensus{i,1} = [specSet{vSets{i,1},1}];
    consensus{i,idxPM} = mean([specSet{vSets{i,1},idxPM}]);

    if isContig==0
        h = hist([specSet{vSets{i,1},4}],[0:3]);
        consensus{i,4} = min(find(h==max(h)))-1;   % -1 because first bin is for charge zero (unknown)
    end;

    % Find the annotation with most occurrences
    maxLen=0; for j=1:size(vSets{i,1},2) if size(specSet{vSets{i,1}(j),idxPep},2)>maxLen  maxLen=size(specSet{vSets{i,1}(j),idxPep},2); end; end;
    if maxLen==0 consensus{i,idxPep}=''; else
        strs = unique(specSet(vSets{i,1},idxPep));     strCounts = zeros(size(strs,1),1);
        for j=1:size(strs,1) strCounts(j)=sum(strcmp(strs{j,1},specSet(vSets{i,1},idxPep))); end;
        consensus{i,idxPep} = strs{ min(find(strCounts==max(strCounts))) ,1};
    end;
end

function [counts, matches, idxPeaks] = aux_ciPeakCoherence(specSet, isContig, windowSize, getIdx, resolution)
% function [counts, matches, idxPeaks] = aux_ciPeakCoherence(specSet, isContig, windowSize, getIdx, resolution)

if nargin<5 resolution=0.1; end;

numSpecs = size(specSet,1);
if isContig maxMass = ceil(max([specSet{:,5}])/resolution); else maxMass = ceil(max([specSet{:,3}])/resolution); end;
window = [-(windowSize/resolution):(windowSize/resolution)];   szWindow = size(window,2);   windowCenter = round(windowSize/resolution)+1;

matches = zeros(maxMass,3);  % First column counts # specs with peaks inside the window, Second/third cols have summed normalized intensities per window (col 2) and per bin (col 3)
if getIdx==1 idxPeaks=cell(maxMass,1); else idxPeaks=[]; end;
for i=1:numSpecs
    if isContig spec = specSet{i,3}; mhMass = specSet{i,5}; else spec = specSet{i,2}; mhMass = specSet{i,3}; end;

    idxKeptPeaks = find(spec(:,1)>windowSize & spec(:,1)<mhMass-windowSize);
    spec = spec(idxKeptPeaks,:);
    if(isempty(spec)) continue; end;

    peaks = round(spec(:,1)/resolution);   peaks = peaks(find(peaks>(windowSize/resolution) & peaks<round((mhMass/resolution)-(windowSize/resolution))));
    peaks = repmat(peaks,1,szWindow)+repmat(window,size(peaks,1),1);

    for j=1:size(peaks,1)
        matches(peaks(j,:),2) = matches(peaks(j,:),2) + spec(j,2) + [-windowSize:resolution:0 -[resolution:resolution:windowSize]]';
        matches(peaks(j,windowCenter),3) = matches(peaks(j,windowCenter),3) + spec(j,2);
        if getIdx==1
            for k=1:size(peaks,2) idxPeaks{peaks(j,k)} = [idxPeaks{peaks(j,k)}; [i idxKeptPeaks(j)]]; end;
        end
    end;

    spec = zeros(maxMass,1);
    spec(peaks)=1;
    matches(:,1) = matches(:,1) + spec;
end;

counts = hist(matches(:,1),[1:numSpecs]);

function newPrms = aux_weightedMerge(prmList, tolerance, scoreType)
% function newPrms = aux_weightedMerge(prmList, tolerance, scoreType)

if isempty(prmList) newPrms = prmList; return; end;
[foo, idxS] = sort(prmList(:,1));
prmsS = prmList(idxS,:);
cur = 1;   newPrms = [];
for j=2:size(prmsS,1)
    if prmsS(j,1) - max(prmsS(cur,1)) <= tolerance
        cur = [cur j];
    else
        sumVal = sum(prmsS(cur,2));  if sumVal==0 sumVal=1; meanMass=mean(prmsS(cur,1)); else meanMass=sum(prmsS(cur,1).*prmsS(cur,2))/sumVal; end;
        if scoreType == 'max'
            newPrms = [newPrms; meanMass max(prmsS(cur,2))];
        else
            newPrms = [newPrms; meanMass sum(prmsS(cur,2))];
        end
        cur = j;
    end
end
sumVal = sum(prmsS(cur,2));  if sumVal==0 sumVal=1; meanMass=mean(prmsS(cur,1)); else meanMass=sum(prmsS(cur,1).*prmsS(cur,2))/sumVal; end;
if scoreType == 'max'
    newPrms = [newPrms; meanMass max(prmsS(cur,2))];
else
    newPrms = [newPrms; meanMass sum(prmsS(cur,2))];
end

function newSpec = aux_mergeConsensusPeaks(spec, lookAhead, resolution)
% function newSpec = aux_mergeConsensusPeaks(spec, lookAhead, resolution)

szSpec = size(spec,1);
peakSets = {};         numSets = 0;
curPeak = spec(1,1);   curSet = [1];
for i=2:szSpec
    if spec(i,1)-curPeak-resolution <= .0001 curSet = [curSet i]; end;
    if i==szSpec | spec(i,1)-curPeak-resolution > .0001
        numSets = numSets+1;
        peakSets = [peakSets; {curSet}];
        curSet = [i];
    end
    curPeak = spec(i,1);
end

w0 = [1:lookAhead];   w1 = round(1/resolution)+w0;   w2 = round(2/resolution)+w0;
newSpec = [];
for i=1:numSets
tmp = spec(peakSets{i},:);
    if size(peakSets{i},2)<=15        scores = scoreWindows(spec(peakSets{i},:),w0,[],[]);
    else if size(peakSets{i},2)<=25   scores = scoreWindows(spec(peakSets{i},:),w0,w1,[]);
        else scores = scoreWindows(spec(peakSets{i},:),w0,w1,w2);
        end
    end
    [intensity, idx] = max(scores);   idx=min(idx);  if intensity<1e-5 continue; end;

    peakIso = idx+10;
    if peakIso <= size(peakSets{i},2)
        if scores(peakIso)==0 % Find nearest peak with score > 0
            all = find(scores>0);   diffs = abs(all-peakIso);  options = find(diffs==min(diffs));
            if max(size(options))==2 options=max(options); end  % Choose farthest peak to cover largest range within tolerance
            peakIso = all(options);
        end
        v = [scores(idx) scores(peakIso)]/sum(scores([idx peakIso]));
        newSpec = [newSpec; spec(peakSets{i}([idx peakIso]'),1) [intensity-scores(peakIso) scores(peakIso)]'];
    else
        newSpec = [newSpec; spec(peakSets{i}(idx),1) intensity];
    end
end

function [scores, scoresIso] = scoreWindows(spec,window1,window2,window3)
% function [scores, scoresIso] = scoreWindows(spec,window1,window2,window3)

szSpec = size(spec,1);   scores = zeros(szSpec,1);    % scoresIso = zeros(szSpec,1);
idx = [window1 window2 window3];   idx = idx(find(idx<=szSpec));
for i=1:szSpec
    if spec(idx(1),2)==0 scores(i)=0; else scores(i) = sum(spec(idx,2)); end % 04/11/01: Require first position != 0
    idx = idx+1;   if idx(size(idx,2))>szSpec idx=idx(1:size(idx,2)-1); end;
end
