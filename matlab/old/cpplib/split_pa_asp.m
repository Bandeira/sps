function [aligns_pa, pairs_asp, aligns_asp, idx_pa] = split_pa_asp(aligns, tolerance)
%
%  Select out alignments where one of the shifts is zero or the parent masses align after the shift (ASP)
%

idx_asp1 = find(abs(aligns(:,3))<=tolerance);   idx_asp2 = find(abs(aligns(:,4))<=tolerance);   idx_asp2 = setdiff(idx_asp2,idx_asp1);
aligns_asp = single([aligns(idx_asp1,4:6); aligns(idx_asp2,[3 5:6])]);
pairs_asp = uint16(aligns([idx_asp1;idx_asp2],1:2));   [pairs_asp,idxS]=sortrows(pairs_asp);   aligns_asp = aligns_asp(idxS,:);

idx_pa = setdiff([1:size(aligns,1)]',[idx_asp1;idx_asp2]);
aligns_pa = aligns(idx_pa,:);
