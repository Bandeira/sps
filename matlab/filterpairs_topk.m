function topk = filterpairs_topk(pairs, k, matchedVerts)

v = unique([pairs(:,1);pairs(:,2)]);
adj = cell(max(v),1);
for i=1:size(pairs,1)
    adj{pairs(i,1)} = [adj{pairs(i,1)}; [pairs(i,5) i]];
    adj{pairs(i,2)} = [adj{pairs(i,2)}; [pairs(i,5) i]];
end

keep = zeros(size(pairs,1),1);
for i=1:size(v,1)
    s = v(i);
    szAdj = size(adj{s},1);
    if szAdj>0
        adj{s} = sortrows(adj{s});
        n = min(k,szAdj);
        idx = adj{s}(szAdj-n+1:szAdj,2);
        keep(idx) = keep(idx) + 1;
    end
end

topk = pairs(keep>=matchedVerts,:);
