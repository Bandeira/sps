function prmLists = mergeConnected(contigs, aligns, vSets, eSets, tolerance);
% function [newContigs, alignsLeftover] = mergeConnected(contigs, aligns, vSets, eSets, tolerance);
%
%  Takes a set of connected components' vertices (vSets) and directed edges (eSets) and
%    merges the contigs accordingly
%
%  prmLists - internal prm matches (1), lower limit prms (2), upper limit prms (3), mean total parent masses (4)
%
%  A PRM appears in the merged contig if it appears in the overlap for _at_least_ one edge,
%    final prm scores are the sum of the scores over all prm-spectra
%  Final filtering step includes a low-pass filter of 1 Da - keep a prm if it's maximal in a 
%    window of 1 Da around it and set m/z = weighted average and score = sum scores
%

numComps = size(vSets,1);
numContigs = size(contigs,1);

% newContigs = cell(numComps, size(contigs,2));   
prmLists = cell(numComps,4);
prmsToKeep = cell(numContigs,1);

for i=1:numComps
    numVerts = size(vSets{i,1},2);   numEdges = size(eSets{i,1},2);
    vertices = vSets{i,1};           edges = aligns(eSets{i,1},:);
    for e=1:numEdges
        [score, prmsList] = scoreOverlapE(contigs{edges(e,1),3}, 0, contigs{edges(e,2),3}, 0, edges(e,3), 6, tolerance);
        prmsToKeep{edges(e,1)} = union(prmsToKeep{edges(e,1)}, prmsList{1,1});
        prmsToKeep{edges(e,2)} = union(prmsToKeep{edges(e,2)}, prmsList{1,2});
    end
    
    % get absolute shifts
    absShifts = getAbsoluteShifts(edges, tolerance);
    
    % Setup the scenario for average prm-scores
    maxPMass = round(10*max(absShifts(vertices)' + [contigs{vertices,5}]));
    prmSpecsCounts = zeros(maxPMass,1);
    for j=1:numVerts
        lims = round(10*(absShifts(vertices(j))+[min(contigs{vertices(j),3}(:,1)) max(contigs{vertices(j),3}(:,1))]));  % Interval where there are prms in spectrum j
        prmSpecsCounts(lims(1):lims(2)) = prmSpecsCounts(lims(1):lims(2)) + ones(lims(2)-lims(1)+1,1);
    end;
        
    prmList = [];  % Compiles complete list of matched prms from all contigs
    for j=1:numVerts
        v = vertices(j);
        prmList = [prmList; weightedMerge([absShifts(v)+contigs{v,3}(prmsToKeep{v},1) contigs{v,3}(prmsToKeep{v},2)], tolerance, 'max')];
    end
    
    % Low pass filter the prmList: should have only correct prefix ions - a different approach is to use scoreOverlap/prmList & rescore around it
%     newPrms = weightedMerge(prmList, tolerance);

    prmLists{i,1} = weightedMerge(prmList, tolerance, 'sum');
    [foo prmIndex] = scoreOverlapE(prmLists{i,1}, 0, prmLists{i,1}, 0, 0, 6, 0);
    prmLists{i,1} = prmLists{i,1}(prmIndex{1,1},:);
    prmLists{i,1}(:,2) = prmLists{i,1}(:,2)./prmSpecsCounts(round(10*prmLists{i}(:,1)));   % average the prm scores
        
%     % Add 'fixed' prms at prm-spectra boundaries (prm-score = max(prm scores)+1), 
%     bounds = [absShifts(vertices)' absShifts(vertices)'+[contigs{vertices,5}]-1 ];

    prmLists{i,2} = absShifts(vertices)';                           prmLists{i,2} = prmLists{i,2}(find(prmLists{i,2}>tolerance));
    
    prmLists{i,3} = absShifts(vertices)'+[contigs{vertices,5}];   
    prmLists{i,4} = mean( prmLists{i,3}( find(prmLists{i,3}>max(prmLists{i,3})-10) ) );
    prmLists{i,3} = prmLists{i,3}(find(prmLists{i,3}<max(prmLists{i,3})-56));

%     % Add the limit prms to the list consensus prms
%     limitsPrms = [prmLists{i,2} prmLists{i,2}+18 prmLists{i,3}-19 prmLists{i,3}-1]';
%     prmLists{i,1} = [prmLists{i,1}; [limitsPrms mean(prmLists{i,1}(:,2))*ones(size(limitsPrms,1),1)]];
%     [foo idx] = sort(prmLists{i,1}(:,1));
%     prmLists{i,1} = prmLists{i,1}(idx,:);
    
%     newPrms = weightedMerge( [ newPrms ; [bounds' max(newPrms(:,2))*ones(size(bounds,2),1)] ] , tolerance);
% 
%     [foo prmIndex] = scoreOverlapE(newPrms, newPrms, 0, 6, 0);
%     newPrms = newPrms(prmIndex{1,1},:);
%     % Enforce prm-spectrum symmetry
end

% newContigs = [newContigs; contigs(setdiff([1:size(contigs,1)],[vSets{:,1}]),:)];

% TODO: relabel contig indices in alignsLeftover


