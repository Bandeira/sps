function [pairs_asp, aligns_asp, ratios_asp, aligns_pa, ratios_pa] = merge_batch_asp_pa(pairs_asp, aligns_asp, ratios_asp, aligns_pa, ratios_pa)
% function [pairs_asp, aligns_asp, ratios_asp, aligns_pa, ratios_pa] = merge_batch_asp_pa(pairs_asp, aligns_asp, ratios_asp, aligns_pa, ratios_pa)
%
%  Removes repeated entries from (pairs/aligns)_asp vs aligns_pa by keeping only the entry with the best summed ratios.
%

% Sort aligns
[foo, idxS] = sortrows(pairs_asp);   pairs_asp = pairs_asp(idxS,:); aligns_asp = aligns_asp(idxS,:); ratios_asp = ratios_asp(idxS,:);
[foo, idxS] = sortrows(aligns_pa(:,1:2));   aligns_pa = aligns_pa(idxS,:); ratios_pa = ratios_pa(idxS,:);

% Remove duplicate entries
idx_asp = 1;   idx_pa = 1;
toKeep_asp = ones(size(pairs_asp,1),1);   toKeep_pa = ones(size(aligns_pa,1),1);
while idx_asp<=size(pairs_asp,1) & idx_pa<=size(aligns_pa,1)
    v1 = double(pairs_asp(idx_asp,1));   v2 = double(pairs_asp(idx_asp,2));
    if v1==aligns_pa(idx_pa,1) & v2==aligns_pa(idx_pa,2)
        if sum(ratios_asp(idx_asp,:))>sum(ratios_pa(idx_pa,:)) toKeep_pa(idx_pa)=0; else toKeep_asp(idx_asp)=0; end;
        idx_pa = idx_pa+1; idx_asp = idx_asp+1; continue;
    end
    
    if v1<aligns_pa(idx_pa,1) | (v1==aligns_pa(idx_pa,1) & v2<aligns_pa(idx_pa,2)) idx_asp = idx_asp+1; else idx_pa = idx_pa+1; end;
end
toKeep_asp = find(toKeep_asp==1);   toKeep_pa = find(toKeep_pa==1);
pairs_asp = pairs_asp(toKeep_asp,:);
aligns_asp = aligns_asp(toKeep_asp,:);
ratios_asp = ratios_asp(toKeep_asp,:);
aligns_pa = aligns_pa(toKeep_pa,:);
ratios_pa = ratios_pa(toKeep_pa,:);
