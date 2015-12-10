function prmLists = mergeConnected3(contigs, prmSpecs, aligns, vSets, eSets, tolerance);
% function prmLists = mergeConnected3(contigs, specSet, aligns, vSets, eSets, tolerance);
%
%  Takes a set of connected components' vertices (vSets) and directed edges (eSets) and
%    merges the contigs _BY_RESCORING_THE_PRMS_ over all the overlapped MS/MS spectra
%
%  prmLists - internal prm matches (1), lower limit prms (2), upper limit prms (3), mean total parent masses (4)
%  prmSpecs - generated using scoreprm5 (with fixed parent mass) for the alpha-synuclein experiments
%
%  A PRM appears in the merged contig if it appears in _any_ overlap,
%    final prm scores are the sum of the scores over all prm-spectra
%  Final filtering step includes a low-pass filter of 1 Da - keep a prm if it's maximal in a 
%    window of 1 Da around it and set m/z = weighted average and score = sum scores
%
%  Similar to mergeConnected but only allows one interpretation per PRM pair. Selection is done
%    in a greedy way by selecting the overlapped prm with highest scores (over all overlapped 
%    PRM spectra) and choosing that interpretation for the pair
%
% 
%  Change log:
%
%    2005/03/02: - Corrected the bug where the scores of the removed PRMs were not removed from each spectrum in prmSpecs (but only from matchScores)
%

numComps = size(vSets,1);   numContigs = size(contigs,1);   
tolRange = round([-tolerance*10:tolerance*10]);

% newContigs = cell(numComps, size(contigs,2));   
prmLists = cell(numComps,4);
prmsToKeep = cell(numContigs,1);

for i=1:numComps
    numVerts = size(vSets{i,1},2);   numEdges = size(eSets{i,1},2);
    vertices = vSets{i,1};           edges = aligns(eSets{i,1},:);

    % get absolute shifts
    absShifts = getAbsoluteShifts(edges, 2*tolerance, 0);
    
    % Data structures to quickly find best matched prm and intervening prms' indices
    maxMass = round(max(absShifts(vertices)'+[contigs{vertices,5}])*10); 
    matchScores = zeros(maxMass,1);   
    for j=1:numVerts   v = vertices(j);

		%    prmSpecs{v,1} - array with the scores for every putative prm (either from scoreprm5 or genPrmSpecs)
		%    prmSpecs{v,2} - minimum prm m/z  
		%    prmSpecs{v,3} - maximum prm m/z
		%    prmSpecs{v,4} - absolute shift
		%    prmSpecs{v,5} - list of prm positions already subtracted from matchScores
		%    prmSpecs{v,6} - minimum PRM score
            
        minPrm = min(find(prmSpecs{v,1}(:,2)>0));    prmSpecs{v,2} = minPrm;
        maxPrm = max(find(prmSpecs{v,1}(:,2)>0));    prmSpecs{v,3} = maxPrm;   
        prmSpecs{v,4} = round(absShifts(v)*10);      prmSpecs{v,5} = [];       prmSpecs{v,6} = min([0;prmSpecs{v,1}(:,2)]);
        matchRange = [prmSpecs{v,4}+minPrm:prmSpecs{v,4}+maxPrm];
        matchScores(matchRange) = matchScores(matchRange) + prmSpecs{v,1}(minPrm:maxPrm,2);  % prmSpecs{v,1} should already be normalized by numSpecs in cluster
    end

    maxIdx = find(matchScores==max(matchScores) & matchScores>0);
prevMaxIdx = -1;  % DEBUG CODE
    while ~isempty(maxIdx) & matchScores(maxIdx(1))>0.1
tmp = [maxIdx/10 matchScores(maxIdx)];   % DEBUG CODE
        % Pick highest scoring sum of matched prms and get the prms' indices in the prm-spectra
        if size(maxIdx,1)>1 
%             maxIdx=weightedMerge([maxIdx/10 matchScores(maxIdx)],tolerance,'max'); maxIdx=round(10*maxIdx(1,1)); 
            maxIdx = selectPRM([maxIdx matchScores(maxIdx)], round(10*tolerance));
        end;
if maxIdx==prevMaxIdx   % DEBUG CODE
    a=1;
end
prevMaxIdx = maxIdx;   % DEBUG CODE
        matchSpecs=0;  
        % Remove resolved prm pairs
        prmScore = 0;
prevMatchScore = matchScores(maxIdx);   % DEBUG CODE
        for j=1:numVerts   v = vertices(j);    prmIdx = maxIdx - prmSpecs{v,4};
%             if prmIdx<=prmSpecs{v,2} | prmIdx>prmSpecs{v,3} continue; end;
            if prmIdx<prmSpecs{v,2} | prmIdx>prmSpecs{v,3} continue; end;
            prmScore   = prmScore + prmSpecs{v,1}(prmIdx,2);   % include prm score if within range -  DEBUG CODE
            if (prmSpecs{v,1}(prmIdx,2)<=0) continue; end;
            matchSpecs = matchSpecs + 1;  % Count in how many specs the prm is matched

            symIdx = round(contigs{v,5}*10-10-prmIdx);
            remIdx = [prmIdx-550:prmIdx+550 symIdx-10:symIdx+10]';
%             remIdx = [prmIdx-10:prmIdx+10 symIdx-10:symIdx+10]';
            
            % Cancel the scores of every position within remIdx that already contributed to matchScores
            remIdx = remIdx(find(remIdx>=prmSpecs{v,2} & remIdx<=prmSpecs{v,3}));
            remIdx = setdiff(remIdx,prmSpecs{v,5});
            
            % PRMs in remIdx can have both positive _and_ negative values - the latter should not be subtracted? Ok, net-effect is <= 0 because of prmSpecs{v,6}
            
            if isempty(remIdx)    % DEBUG CODE
                continue; 
            end;
            prmSpecs{v,5} = [prmSpecs{v,5}; remIdx];
            matchScores(prmSpecs{v,4}+remIdx) = matchScores(prmSpecs{v,4}+remIdx) - prmSpecs{v,1}(remIdx,2);

            prmSpecs{v,1}(remIdx,2) = prmSpecs{v,6};  % Set scores as if prmSpecs{v,1} never had a prm at this position
            matchScores(prmSpecs{v,4}+remIdx) = matchScores(prmSpecs{v,4}+remIdx) + prmSpecs{v,6};
        end
        if abs(prmScore-prevMatchScore)>1e-5    % DEBUG CODE
            fprintf(1,'prmScore~=prevMatchScore!!  prmScore == %.1f, prevMatchScore == %.1f\n',prmScore,prevMatchScore);
        end
        if matchSpecs>1 prmLists{i,1} = [prmLists{i,1}; maxIdx/10 prmScore]; end;
        maxIdx = find(matchScores==max(matchScores) & matchScores>0);
        tmp = prmLists{i,1};
    end;

%     % Setup the scenario for average prm-scores
%     maxPMass = round(10*max(absShifts(vertices)' + [contigs{vertices,5}]));
%     prmSpecsCounts = zeros(maxPMass,1);
%     for j=1:numVerts
%         lims = round(10*(absShifts(vertices(j))+[min(contigs{vertices(j),3}(:,1)) max(contigs{vertices(j),3}(:,1))]));  % Interval where there are prms in spectrum j
%         prmSpecsCounts(lims(1):lims(2)) = prmSpecsCounts(lims(1):lims(2)) + ones(lims(2)-lims(1)+1,1);
%     end;
%         
%     prmLists{i,1} = weightedMerge(prmLists{i,1}, tolerance, 'sum');
%     [foo prmIndex] = scoreOverlapE(prmLists{i,1}, 0, prmLists{i,1}, 0, 0, 6, 0);
%     prmLists{i,1} = prmLists{i,1}(prmIndex{1,1},:);
%     prmLists{i,1}(:,2) = prmLists{i,1}(:,2)./prmSpecsCounts(round(10*prmLists{i}(:,1)));   % average the prm scores

    prmLists{i,2} = absShifts(vertices)';                           prmLists{i,2} = prmLists{i,2}(find(prmLists{i,2}>tolerance));
    
    prmLists{i,3} = absShifts(vertices)'+[contigs{vertices,5}];   
    prmLists{i,4} = mean( prmLists{i,3}( find(prmLists{i,3}>max(prmLists{i,3})-10) ) );
    prmLists{i,3} = prmLists{i,3}(find(prmLists{i,3}<max(prmLists{i,3})-56));

end


function prmMass = selectPRM(prmSpec, tolerance)
%
%  Selects the center-most high-scoring PRM from the first group of PRMs. A group of PRMs is a sequence
%    of PRMs where the distance between any two consecutive PRM masses is <= tolerance
%
%    The selected PRM is _guaranteed_ to be one of the PRMs in prmSpec.
%
%  NOTE: tolerance and mass values in prmSpec should be integer values.
%

% Find first group
numPRMs = size(prmSpec,1); groupEnd=1;
while groupEnd<numPRMs & prmSpec(groupEnd+1,1)-prmSpec(groupEnd,1)-tolerance<=1e-5  groupEnd=groupEnd+1; end;
prmSpec = prmSpec(1:groupEnd,:);
numPRMs = groupEnd;
if numPRMs==1 prmMass=prmSpec(1,1); return; end

% Low pass filter
arraySpec = zeros(prmSpec(numPRMs,1)-prmSpec(1,1)+2*tolerance+1,1); % Array version of prmSpec for low pass filter
arrayIndices = prmSpec(:,1)-prmSpec(1,1)+tolerance+1;               % Indices in arraySpec for every PRM in prmSpec
arraySpec(arrayIndices) = prmSpec(:,2);
lpSpec = prmSpec;
for i=1:numPRMs lpSpec(i,2)=sum(arraySpec(arrayIndices(i)-tolerance:arrayIndices(i)+tolerance)); end;
weightedCenter = sum( lpSpec(:,1).*(lpSpec(:,2)/sum(lpSpec(:,2))) );
prmSpec = lpSpec(find(lpSpec(:,2)==max(lpSpec(:,2))),:);
if size(prmSpec,1)==1 prmMass=prmSpec(1,1); return; end

% Find center-most PRM
distToCenter = abs(prmSpec(:,1)-weightedCenter);
prmSpec = prmSpec(find(distToCenter==min(distToCenter)),:);
prmMass = prmSpec(1,1); 
