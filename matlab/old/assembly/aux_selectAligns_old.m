function subset = aux_selectAligns(aligns, pepsToKeep, pepsToGo, idxOffset)
%
%   Auxiliary function for assembleC and assemble1. Selects and renumbers the alignments 
%     kept from one assembly iteration to the next.
%

    szPepsToKeep = size(pepsToKeep,2);
    
    newIdx = zeros(max([pepsToKeep pepsToGo]), 1);
    newIdx(pepsToKeep) = idxOffset + [1:szPepsToKeep];  % renumber the alignments
    
    subset = zeros(szPepsToKeep*(szPepsToKeep-1), size(aligns,2));
    cur = 1;
    for i=1:size(aligns,1)
        if newIdx(aligns(i,1))>0 & newIdx(aligns(i,2))>0
            subset(cur,:) = [newIdx(aligns(i,1)) newIdx(aligns(i,2)) aligns(i,3:5)];
            cur=cur+1;
        end
    end
