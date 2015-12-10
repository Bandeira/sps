function pos = getAssembledRange(tagsMatch, pMatch, seqsMatched, protIdx, components)
% function pos = getAssembledRange(tagsMatch, pMatch, seqsMatched, protIdx, components)
%
%  Determines the ranges covered by assembled peptides, according to the individual spectrum IDs from database search.
%   The purpose is to visualize the assembled regions and not just the sequenced regions (as obtained from filterLongestTag+findTags).
%
%  tagsMatch/pMatch - Location of all individual peptides on the target proteins (as returned by findTags)
%  seqsMatched - indices and positions of the contigs whose assembly coverage ranges are to be determined (as returned by findTags of the recovered seqs, pMatch{i,4})
%  protIdx - target protein index as in tagsMatch
%

% if nargin<=4 components = load_binListArray('../assembly/components.bla', 'int'); for c=1:size(components,1) if ~isempty(components{c}) components{c}=components{c}+1; end; end; end; % +1 because components have zero-based spectrum indices
components = components(seqsMatched(:,1));  numComps = size(components,1);   
pos = seqsMatched(:,2:3);   % Initialize with the range where the de-novo sequence matched

% matchedProteins = [];
for c=1:numComps
    specIdx = components{c}';
%     pIdx = tagsMatch(specIdx);   
%     pIdx=pIdx(find(pIdx>0));   pIdxU = unique(pIdx);   if isempty(pIdxU) continue; end;
%     if size(pIdxU,1)>1 h=hist(pIdx,pIdxU); i=min(find(h==max(h))); p=pIdxU(i); else p=pIdxU; end;
%     matchedProteins = [matchedProteins p];
%     specIdx = specIdx( find(tagsMatch(specIdx)==p) );  % Consider only peptides that matched the most-matched protein (p)
%     pIdx = find([pMatch{:,1}]==p);

    specIdx = specIdx( find(tagsMatch(specIdx)==protIdx) );  if isempty(specIdx) continue; end;  % Consider only peptides that matched the target protein (protIdx)
    pIdx = find([pMatch{:,1}]==protIdx);   
    specPos = pMatch{pIdx,4}( csps_matchCols(pMatch{pIdx,4}(:,1),specIdx) ,2:3);   

    covered = zeros(max(specPos(:,2)),1);  szCovered = size(covered,1);   % Find the largest consecutive interval covered by assembled spectra
    for s=1:size(specPos,1) covered(specPos(s,1):specPos(s,2))=1; end;
%     covered = cumsum(covered);    szCovered = size(covered,1);    covered(find(covered(1:szCovered-1)==covered(2:szCovered))+1)=0;
    idxS = find(covered(1:szCovered-1)==0 & covered(2:szCovered)==1)+1;   if covered(1)==1 idxS=[1; idxS]; end;
    idxE = find(covered(1:szCovered-1)==1 & covered(2:szCovered)==0);     if covered(szCovered)==1 idxE=[idxE; szCovered]; end;
    for i=1:size(idxS,1) covered(idxS(i):idxE(i)) = cumsum(covered(idxS(i):idxE(i))); end;
    endPos = min(find(covered==max(covered)));   startPos = max(find(covered(1:endPos)==1));
    
    if startPos>pos(c,2) | endPos<pos(c,1)  % Intervals don't overlap, keep the largest
        if endPos-startPos > pos(c,2)-pos(c,1) pos(c,:)=[startPos endPos]; end;
    else
        pos(c,:) = [min(pos(c,1),startPos) max(pos(c,2),endPos)];
    end
end;

% matchedProteins = unique(matchedProteins);   fprintf(1,'Matched %d different proteins: ',size(matchedProteins,2));
% for p=1:size(matchedProteins,2) fprintf(1,'%d ',matchedProteins(p)); end; fprintf(1,'\n');