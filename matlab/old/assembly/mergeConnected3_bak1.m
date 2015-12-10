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

numComps = size(vSets,1);   numContigs = size(contigs,1);   
tolRange = round([-tolerance*10:tolerance*10]);

% newContigs = cell(numComps, size(contigs,2));   
prmLists = cell(numComps,4);
prmsToKeep = cell(numContigs,1);

for i=1:numComps
    numVerts = size(vSets{i,1},2);   numEdges = size(eSets{i,1},2);
    vertices = vSets{i,1};           edges = aligns(eSets{i,1},:);

    % get absolute shifts
    absShifts = getAbsoluteShifts(edges, 2*tolerance);
    
    % Data structures to quickly find best matched prm and intervening prms' indices
    maxMass = round(max(absShifts(vertices)'+[contigs{vertices,5}])*10); 
    matchScores = zeros(maxMass,1);   
    for j=1:numVerts   v = vertices(j);

		%    prmSpecs{v,1} - array with the scores for every putative prm (either from scoreprm5 or genPrmSpecs)
		%    prmSpecs{v,2} - minimum prm m/z  
		%    prmSpecs{v,3} - maximum prm m/z
		%    prmSpecs{v,4} - absolute shift
		%    prmSpecs{v,5} - list of prm positions already subtracted from matchScores
            
        minPrm = min(find(prmSpecs{v,1}(:,2)>0));    prmSpecs{v,2} = minPrm;
        maxPrm = max(find(prmSpecs{v,1}(:,2)>0));    prmSpecs{v,3} = maxPrm;   
        prmSpecs{v,4} = round(absShifts(v)*10);      prmSpecs{v,5} = [];
        matchRange = [prmSpecs{v,4}+minPrm:prmSpecs{v,4}+maxPrm];
        matchScores(matchRange) = matchScores(matchRange) + prmSpecs{v,1}(minPrm:maxPrm,2);  % prmSpecs{v,1} should already be normalized by numSpecs in cluster
    end

    maxIdx = find(matchScores==max(matchScores) & matchScores>0);
    while ~isempty(maxIdx) & matchScores(maxIdx(1))>0.1
        % Pick highest scoring sum of matched prms and get the prms' indices in the prm-spectra
        if size(maxIdx,1)>1 
            maxIdx=weightedMerge([maxIdx/10 matchScores(maxIdx)],tolerance,'max'); maxIdx=round(10*maxIdx(1,1)); 
        end;
        matchSpecs=0;  
        % Remove resolved prm pairs
        prmScore = 0;
        for j=1:numVerts   v = vertices(j);    prmIdx = maxIdx - prmSpecs{v,4};
%             if prmIdx<=prmSpecs{v,2} | prmIdx>prmSpecs{v,3} continue; end;
            if prmIdx<prmSpecs{v,2} | prmIdx>prmSpecs{v,3} continue; end;
            prmScore   = prmScore + prmSpecs{v,1}(prmIdx,2);   % include prm score if within range
            if (prmSpecs{v,1}(prmIdx,2)<=0) continue; end;
            matchSpecs = matchSpecs + 1;  % Count in how many specs the prm is matched

            symIdx = round(contigs{v,5}*10-10-prmIdx);
            remIdx = [prmIdx-550:prmIdx+550 symIdx-10:symIdx+10]';
%             remIdx = [prmIdx-10:prmIdx+10 symIdx-10:symIdx+10]';
            
            % Cancel the scores of every position within remIdx that already contributed to matchScores
            remIdx = remIdx(find(remIdx>=prmSpecs{v,2} & remIdx<=prmSpecs{v,3}));
            remIdx = setdiff(remIdx,prmSpecs{v,5});
            if isempty(remIdx) continue; end;
            prmSpecs{v,5} = [prmSpecs{v,5}; remIdx];
            matchScores(prmSpecs{v,4}+remIdx) = matchScores(prmSpecs{v,4}+remIdx) - prmSpecs{v,1}(remIdx,2);
            % and then add negative scores as if prmSpecs{v,1} never had a prm at this position - prmSpecs{v,1}(1,2) should have lowest score
%             matchScores(prmSpecs{v,4}+remIdx) = matchScores(prmSpecs{v,4}+remIdx) + prmSpecs{v,1}(1,2);
            matchScores(prmSpecs{v,4}+remIdx) = matchScores(prmSpecs{v,4}+remIdx) + min(prmSpecs{v,1}(:,2));
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

