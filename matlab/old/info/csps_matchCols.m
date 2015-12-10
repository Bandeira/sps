function idx=matchCols(c1,c2)
% function idx=matchCols(c1,c2)
%
%  idx is the index in c1 of elements in c2
%

matches = zeros(size(c1,1),1);
for i=1:size(c2,1)
    matches(find(c1==c2(i)))=1;
end
idx=find(matches==1);
