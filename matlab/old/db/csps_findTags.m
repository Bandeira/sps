function [tagsMatch, protsMatch, protsGroups] = findTags(tags, db)
% function [tagsMatch, protsMatch, protsGroups] = findTags(tags, db)
%
%  Finds which proteins are matched by which tags, reranks proteins by decreasing number of matched tags
%   and finally reassigns tag matches by protein ranks (highest rank first).
%
%  tags - cell columns of peptide strings (e.g. as returned by filterLongestTag, col.2)
%  db   - cell columns of protein strings (e.g. as returned by fasta_getentries, col.2)
%
%  tagsMatch  - column vector, i-th position has protein index for i-th tag (col.1), start/end match pos (cols 2/3)
%  protsMatch - cell array, matched protein indexes sorted by decreasing rank.
%                 Protein index (col.1), Number of matched tags (col.2), Protein length (col.3), List of matched tags (col 4.1) and their positions (cols 4.2-4.3), sorted by increasing tag index (col.4)
%  protsGroups - cell array with groups of matched proteins. Each protsGroups(i,:) contains a set of protein indices (protsGroups{i,1}(j)) and list of matched peptides (protsGroups{i,2}{j}).
%                 Formally, a group is a connected component in the bipartite peptide/protein graph.
%

szTags = size(tags,1);         szDb = size(db,1);
tagsMatch = zeros(szTags,3);   protsMatch = {};
protsMatched = cell(szTags,1);   % Keeps track of what proteins each tag matched
tagsMatched = cell(szDb,1);      % Keeps track of what tags matched the current protein

% Find tags over protein sequences
numPairs = 0;
for tagIdx=1:szTags
    for protIdx=1:szDb
        idx = strfind(db{protIdx},tags{tagIdx});
        if ~isempty(idx)
            protsMatched{tagIdx} = [protsMatched{tagIdx} protIdx];
            tagsMatched{protIdx} = [tagsMatched{protIdx} tagIdx];
            numPairs = numPairs+1;
        end;
    end
end

% Build protein groups
pepProtPairs = zeros(numPairs,2);   numPairs=0;
for tagIdx=1:szTags
    numPeps = length(protsMatched{tagIdx});
    pepProtPairs(numPairs+1:numPairs+numPeps,1) = tagIdx;
    pepProtPairs(numPairs+1:numPairs+numPeps,2) = szTags+protsMatched{tagIdx}';   % +szTags guarantees that peptide/protein indices are different
    numPairs = numPairs+numPeps;
end;
[vSets,eSets] = csps_assembleE(pepProtPairs);
numSets = size(eSets,1);    protsGroups = cell(numSets,2);
for setIdx=1:numSets
    curPairs = pepProtPairs(eSets{setIdx},:);
    protsGroups{setIdx,1} = unique(curPairs(:,2));    numProts = length(protsGroups{setIdx,1});
    protsGroups{setIdx,2} = cell(numProts,1);
    for pIdx=1:numProts
        idx = find(curPairs(:,2)==protsGroups{setIdx,1}(pIdx));
        protsGroups{setIdx,2}{pIdx} = curPairs(idx,1);
    end
    protsGroups{setIdx,1} = protsGroups{setIdx,1} - szTags;   % Correct the protein indices
end

% Find proteins ranks and assign tag matches depending on protein ranks
h=hist([protsMatched{:}],[1:szDb]);   processedTags = [];
while max(h)>0
    highestRankIdx = min(find(h==max(h)));    protsMatch = [protsMatch; {highestRankIdx, h(highestRankIdx), size(db{highestRankIdx},2), []}];
    curMatches = setdiff(tagsMatched{highestRankIdx},processedTags);    protsMatch{size(protsMatch,1),4} = curMatches';
    processedTags = [processedTags curMatches];
    tagsMatch(curMatches,1) = highestRankIdx;
    for tagIdx=1:size(curMatches,2)
        protsMatched{curMatches(tagIdx)} = [];  % Block matched tags from voting for the next highest rank protein
    end;
    h=hist([protsMatched{:}],[1:szDb]);
end

% Get tag positions on the matched protein sequence
for protIdx=1:size(protsMatch,1)
    pos = zeros(protsMatch{protIdx,2},2);
    for tagIdx=1:protsMatch{protIdx,2}
        tagPos = min(strfind(db{protsMatch{protIdx,1}},tags{protsMatch{protIdx,4}(tagIdx)}));
        pos(tagIdx,:) = [tagPos tagPos+size(tags{protsMatch{protIdx,4}(tagIdx)},2)-1];
    end
    tagsMatch(protsMatch{protIdx,4},2:3) = pos;
    protsMatch{protIdx,4} = [protsMatch{protIdx,4} pos];
end