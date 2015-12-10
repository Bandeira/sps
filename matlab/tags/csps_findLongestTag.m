function tagLen = findLongestTag(tag, peakTol)
% function tagLen = findLongestTag(tag, peakTol)
%
% Find the longest tag in a de-novo interpretation allowing for one jump of size 2
%

global AAmasses AAletters;
peakTolInt = round(10*peakTol);   tolRange = [-peakTolInt:peakTolInt];   szTolRange = size(tolRange,2);
jumps1 = unique(repmat(round(10*AAmasses),1,szTolRange)+repmat(tolRange,size(AAmasses,1),1));   jumps1ok = zeros(max(jumps1),1);   jumps1ok(jumps1)=1;   clear jumps1;
jumps2 = csps_getjumps(AAmasses,AAletters,2);   jumps2 = unique(repmat(round(10*jumps2),1,szTolRange)+repmat(tolRange,size(jumps2,1),1));   jumps2ok = zeros(max(jumps2),1);   jumps2ok(jumps2)=1;   clear jumps2;

masses = round(10*getmasses3(tag,'',[]));   szMasses = size(masses,2);   massType = zeros(1,szMasses);
for i=1:szMasses
    if masses(i)<=size(jumps1ok,1) & jumps1ok(masses(i))==1 massType(i)=1; continue; end;
    if masses(i)<=size(jumps2ok,1) & jumps2ok(masses(i))==1 massType(i)=2; end;
end;

% Find longest tag, allowing for 1 jump of size 2
tagLen=0;
tagStart = min(find(massType==1));  if isempty(tagStart) if ~isempty(find(massType==2)) tagLen=2; end; return; end;
tagEnd = tagStart;   jumpPos = 0;
idxJumps = max(find(find(massType==2)<tagStart));  if ~isempty(idxJumps) tagStart=idxJumps; jumpPos=tagStart; end; % Allow the tag to start at the last jump of size 2 before the first amino acid
while tagEnd<=szMasses
    if massType(tagEnd)==2 
        if jumpPos>0 tagLen=max(tagLen,tagEnd-tagStart+1); tagStart = jumpPos+1; end;  % +1 in tagLen is for jump of size 2 included in the tag
        jumpPos=tagEnd;
    end
    if massType(tagEnd)==0
        if jumpPos>0 tagLen=max(tagLen,tagEnd-tagStart+1); else tagLen=max(tagLen,tagEnd-tagStart); end;  % +1 in tagLen is for jump of size 2 included in the tag
        tagStart = tagEnd+1;   jumpPos=0;
    end
    tagEnd=tagEnd+1;
end
if jumpPos>0 tagLen=max(tagLen,tagEnd-tagStart+1); else tagLen=max(tagLen,tagEnd-tagStart); end;  % +1 in tagLen is for jump of size 2 included in the tag
