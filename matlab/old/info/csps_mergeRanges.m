function pos = mergeRanges(pos)
% function pos = mergeRanges(pos)
%
% Merges overlapping or immediately adjacent intervals.
%
% pos - 2-col vector of [start,end] positions
%

maxPos = max(pos(:,2));
covered = zeros(maxPos,1);
for i=1:size(pos,1)
    covered(pos(i,1):pos(i,2))=1;
end;

if covered(1)==1 startPos=1; else startPos=[]; end;
startPos = [startPos; find(covered(1:maxPos-1)==0 & covered(2:maxPos)==1)+1];
endPos = [find(covered(1:maxPos-1)==1 & covered(2:maxPos)==0); maxPos];  % Note that maxPos is guaranteed to be the end of the rightmost interval

pos = [startPos endPos];
