function alignsOut = filterDirection(contigs, aligns, tolerance)
%
%  Purpose: merge connected components in the alignment graph based on changing the prm-scores
%
%  Problem: Sometimes leaves a zig-zag in the alignment graph (both correct and symetric shifts)
%
%  NOTE: OLDER, OUTDATED and PROBLEMATIC version - use filterDirection2

return;

% 1 - Given alignments, get connected components
[vSets, eSets] = assembleE(aligns);
numComps = size(vSets,1);

edgesPerVertex = cell(size(contigs,1),1);
for i=1:size(aligns,1)
    edgesPerVertex{aligns(i,1),1} = [edgesPerVertex{aligns(i,1),1} i];
    edgesPerVertex{aligns(i,2),1} = [edgesPerVertex{aligns(i,2),1} i];
end;

keptEdges = cell(numComps,1);
alignsOut = [];
for i=1:numComps
    curAligns = aligns(eSets{i,1},:);
    
    bestEdgeIdx = find( curAligns(:,4)==max(curAligns(:,4)) );  
    bestEdgeIdx = bestEdgeIdx(1);            
    bestEdge = curAligns(bestEdgeIdx,:);
    
    v1 = bestEdge(1,1);   v2 = bestEdge(1,2);
    
    [score, prmIndex] = scoreOverlapE(contigs{v1,3}, contigs{v2,3}, bestEdge(1,3), 6, tolerance);
    
%     contigs{v1,3}(prmIndex{1,1},2) = contigs{v1,3}(prmIndex{1,1},2) + contigs{v2,3}(prmIndex{1,2},2);
%     contigs{v2,3}(prmIndex{1,2},2) = contigs{v2,3}(prmIndex{1,2},2) + contigs{v1,3}(prmIndex{1,1},2);

    [contigs(v1,:), contigs(v2,:)] = rescoreContigs(contigs(v1,:), contigs(v2,:), bestEdge(1,3), prmIndex);

    % Remove other edges between v1 & v2
    keptEdges{i,1} = bestEdgeIdx;
    idx = find((curAligns(:,1)==v1 & curAligns(:,2)==v2) | (curAligns(:,1)==v2 & curAligns(:,2)==v1));
    curAligns(idx,:) = zeros(max(size(idx)),5);
    
    processed = [v1 v2];   unprocessed = setdiff(vSets{i,1}, [v1 v2]);
    
    while size(unprocessed,2)>0
        nextVertices = setdiff(unique(aligns([edgesPerVertex{processed,1}],[1 2])), processed);
        v = nextVertices(1);

        vToProc = intersect(unique(aligns([edgesPerVertex{v,1}],[1 2])), processed);
        for j=1:size(vToProc,2)
            edgesIdx = find((curAligns(:,1)==v & curAligns(:,2)==vToProc(j)) | (curAligns(:,1)==vToProc(j) & curAligns(:,2)==v));
            edges = curAligns(edgesIdx,:);
            numEdges = size(edges,1);
            curAligns(edgesIdx,:) = zeros(numEdges,5);
            
            % Rescore the edges
            scores = zeros(numEdges,1);   prmLists = cell(numEdges,1);
            for e=1:numEdges
                [scores(e) prmLists{e,1}] = scoreOverlapE(contigs{edges(e,1),3}, contigs{edges(e,2),3}, edges(e,3), 6, tolerance);
            end;
            
            bestIdx = find(scores==max(scores));
            if max(size(bestIdx))>1 & sum(edges(bestIdx,3))>0
                fprintf(1,'Hmmm.... Better look here. \n');
            end
            bestIdx = bestIdx(1);   
            keptEdges{i,1} = [keptEdges{i,1} edgesIdx(bestIdx)];
            bestEdge = edges(bestIdx,:);
            
            [contigs(bestEdge(1,1),:) contigs(bestEdge(1,2),:)] = rescoreContigs(contigs(bestEdge(1,1),:),contigs(bestEdge(1,2),:), bestEdge(1,3), prmLists{bestIdx,1});
        end
        
        processed = [processed v];
        unprocessed = setdiff(unprocessed,[v]);
    end       
    % All spectra are re-scored. How to define the consensus prm-spectrum?
    
    alignsOut = [alignsOut; aligns(eSets{i,1}(keptEdges{i,1}),:)];
end


function [contig1, contig2] = rescoreContigs(contig1, contig2, shift, prmLists)
    contig1{1,3}(prmLists{1,1},2) = contig1{1,3}(prmLists{1,1},2) + contig2{1,3}(prmLists{1,2},2);
    contig2{1,3}(prmLists{1,2},2) = contig1{1,3}(prmLists{1,1},2) + contig2{1,3}(prmLists{1,2},2);

    % ADD SCORE FOR THE SPECTRUM EXTREMES -> 2*max score AFTER the sum
    if shift<0 fprintf(1,'Shift should always be >= 0!!!!!\n'); end;

    if size(contig1{1,3},2)==3 extra=0; else extra=[]; end
    contig1{1,3} = [contig1{1,3}; [shift 2*max(contig1{1,3}(:,2)) extra]];
    if shift+contig2{1,5} < contig1{1,5}  contig1{1,3} = [contig1{1,3}; [shift+contig2{1,5} 2*max(contig1{1,3}(:,2)) extra]];  end;
    [foo idx] = sort(contig1{1,3}(:,1));   contig1{1,3} = contig1{1,3}(idx,:);

    if contig1{1,5}-shift < contig2{1,5}  
        contig2{1,3} = [contig2{1,3}; [contig1{1,5}-shift 2*max(contig2{1,3}(:,2)) extra]];  
        [foo idx] = sort(contig2{1,3}(:,1));   contig2{1,3} = contig2{1,3}(idx,:);
    end;
    
    % Remove mirror images ?
    
