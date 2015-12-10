function [tags,counts] = findLongestTags2(specSet, szMaxJump, expInt, peakTol, mods, modMasses, addModMass)
% function [tags,counts] = findLongestTags2(specSet, szMaxJump, expInt, peakTol, mods, modMasses, addModMass)
%
%  Locates and extracts the longest possible tag for each spectrum in the dataset allowing _one_ jump of size up to szMaxJump
%
%  szMaxJump - number of amino acids allowed in a jump
%
%  tags(i,:) - size of the longest tag (col.1), number of jumps in the tag (col.2)
%    or
%  tags(i,:)   - max tag using one jump of [0:szMaxJump] amino-acids
%  counts - number of valid locations for one jump of [0:szMaxJump] amino-acids. Equal to # gapped seqs (1 gap) to generate
%
%
%  NOTE: Assumes that spectra are PRM spectra
%

if size(specSet,2)==5 idxSpec=2; idxPep=5; else idxSpec=3; idxPep=7; end;

canMiss = szMaxJump-1;  % Number of PRMs that can be missed

numSpecs = size(specSet,1);   tags = zeros(numSpecs,szMaxJump);  counts = zeros(numSpecs,1);
for i=1:numSpecs
    szSpec = size(specSet{i,idxSpec},1);   if szSpec<2 || isempty(specSet{i,idxPep}) continue; end;
    
    masses = getmasses4(specSet{i,idxPep}, mods, modMasses, addModMass);  numMasses = size(masses,2);
    if isempty(expInt) | expInt(i,1)>=expInt(i,2) masses = cumsum(masses); else masses=[cumsum(masses(numMasses:-1:1))+18]; end;
    
    diffs = abs(repmat(specSet{i,idxSpec}(:,1),1,numMasses) - repmat(masses,szSpec,1));
    
    % Convert diffs to set of matched pairs of peaks
    allMatches = find(diffs<=peakTol);   if isempty(allMatches) continue; end;
    rows = mod(allMatches,szSpec);       rows(find(rows==0))=szSpec;        cols = ceil(allMatches/szSpec);
    
    pairs = [rows cols];   if isempty(find(pairs(:,2))==numMasses) pairs = [pairs; szSpec+1 numMasses]; end;  % Parent mass is always matched
    tagLen = zeros(2,szMaxJump);  % Line 1 has cur tag length and line 2 has maximal tag length
    szAllowedJumps = [1:szMaxJump];    premium = ones(1,szMaxJump);
    for j=1:size(pairs,1)
        if pairs(j,1)==1 & pairs(j,2)==1
            tagLen(1,:)=2*ones(1,szMaxJump);   % if both first peaks match we collect a bonus for the 0-th peak
        else
            if j==1 
                tagLen(1,:)=premium;
                if pairs(j,1)~=1 continue; end;  % Spectrum already has an error peak before this pair. Restart here.
                tagLen(1,:)=tagLen(1,:)+(pairs(j,2)).*(szAllowedJumps>=pairs(j,2));  % Jumped somewhere into the middle of the spectrum
                counts(i) = counts(i)+1;
            else
                if pairs(j,1)==pairs(j-1,1)+1 & pairs(j,2)==pairs(j-1,2)+1 tagLen(1,:)=tagLen(1,:)+premium; continue; end;  % consecutive true peaks in both
                if pairs(j,1)>pairs(j-1,1)+1 tagLen(2,:)=max(tagLen); tagLen(1,:)=premium; counts(i) = counts(i)+1; continue; end;   % Error in the spectrum
                test = szAllowedJumps>=(pairs(j,2)-pairs(j-1,2));  % Jump
                copyFromLower = find(test(2:szMaxJump)==1);   reset = find(test==0);
                tagLen(2,:)=max(tagLen);
                tagLen(1,1+copyFromLower) = tagLen(1,1);   tagLen(1,reset)=ones(size(reset));
                tagLen(1,:) = tagLen(1,:) + (pairs(j,2)-pairs(j-1,2))*test;
                counts(i) = counts(i)+1;
            end
        end
    end

    tags(i,:) = max(tagLen);
end
