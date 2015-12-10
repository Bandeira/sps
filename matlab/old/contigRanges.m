function contigRanges(data, contigs, starts)

ranges = cell(size(starts,2),1);
begins = zeros(size(starts,2),1);

for c=1:size(starts,2)
    s = starts(c);
    if size(find(contigs==s),1)>0 continue; end;
    range = [s];
    posS=data{starts(c),1};
    begins(c)=posS;
    peps=0;
    while contigs(s)>0
        s = contigs(s);
        range = [range s];
        peps=peps+1;
    end;
    posE = data{s,1} + size(data{s,2},2) - 1;
%    ranges(c,:) = [posS peps posE];
    ranges{c} = range;
end;

idx = find(begins>0);
ranges = ranges(idx);
begins = begins(idx);

[foo idx] = sort(begins);
ranges = ranges(idx,:);

for c=1:size(ranges,1)
    range = ranges{c};
    for r=1:size(range,2)
        fprintf(1,'[%d](%d,%d) - ',range(r), data{range(r),1}, data{range(r),1}+size(data{range(r),2},2)-1);
    end;
    fprintf(1,'\n');
%    fprintf(1,'%d-(%d)-%d\n',ranges(c,:));
end;