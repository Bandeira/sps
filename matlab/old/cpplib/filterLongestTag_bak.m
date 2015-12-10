function tags = filterLongestTag(report, minTagSize, maxJump, maxMod)
% function tags = filterLongestTag(report, minTagSize, maxJump, maxMod)
%
%  Find spectra whose annotation has a tag of at least minTagSize after ignoring abs(mods)<=maxMod Da and jumps<=maxJump.
%    report is as returned by matchma.
%
%  tags(:,1) - Spectrum index
%  tags(:,2) - Longest tag matching selection criteria
%  tags(:,3) - Length of longest tag
%  tags(:,4) - Total mass before the tag
%  tags(:,5) - Total mass after the tag
%

tags = {};   numAnnots = size(report,1);
for i=1:numAnnots
    tp = [report{i}{:,9}];   idxBest = min(find(tp==max(tp)));   if isempty(report{i}{idxBest,8}) continue; end;
    [curTag, massPrefix, massSuffix] = lcl_findLongestTag(report{i}{idxBest,8}, minTagSize, maxJump, maxMod);
    if size(curTag,2)>=minTagSize tags = [tags; {i, curTag, size(curTag,2), massPrefix, massSuffix}]; end;
end

function tag = lcl_findLongestTag(annotation, minTagSize, maxJump, maxMod)

tag = annotation;

% Ignore start/end jumps
idxS=find(tag=='{');   idxE=find(tag=='}');
for i=1:size(idxS,2) tag(idxS(i):idxE(i))='.'; end;

% Ignore jumps<=maxJump
idxS=find(tag=='(');   idxE=find(tag==')');
for i=1:size(idxS,2)
    if idxE(i)-idxS(i)<=maxJump+1 tag([idxS(i) idxE(i)])=' ';
    else tag(idxS(i):idxE(i))='.'; end;
end

% Remove mods
idxS=find(tag=='[');   idxE=find(tag==']');
for i=1:size(idxS,2)
    strMod = tag(idxS(i)+1:idxE(i)-1);   sepIdx = find(strMod==',');   if isempty(sepIdx) tag(idxS(i):idxE(i))='.'; continue; end;
    modMass = str2num(strMod(sepIdx(2)+1:size(strMod,2)));
    if abs(modMass)<=maxMod & sepIdx(2)-sepIdx(1)<=minTagSize tag([idxS(i):idxS(i)+sepIdx(1) idxS(i)+sepIdx(2):idxE(i)])=' ';
    else tag(idxS(i):idxE(i))='.'; end;
end

% Remove spaces
tag = tag(find(tag~=' '));   szTag = size(tag,2);
idxS = find(tag(1:szTag-1)=='.' & tag(2:szTag)~='.');  idxS=idxS+1;  if tag(1)~='.' idxS=[1 idxS]; end;  % All tag-start positions
idxE = find(tag(2:szTag)=='.' & tag(1:szTag-1)~='.');  if tag(szTag)~='.' idxE=[idxE szTag]; end;        % All tag-end positions

if isempty(idxS) tag=''; else
    bestIdx = min(find((idxE-idxS)==max(idxE-idxS)));
    tag = tag(idxS(bestIdx):idxE(bestIdx));
end;
