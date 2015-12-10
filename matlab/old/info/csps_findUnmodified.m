function cPeps = findUnmodified(peptides)
% function cPeps = findUnmodified(peptides)
%
%  Generates the unmodified version of every modified peptide and locates the index of
%    the unmodified peptide in peptides (if present)
%
%  cPeps{i,1} - unmodified peptide (or same peptide if unmodified)
%  cPeps{i,2} - number of detected modifications
%

numPeps = size(peptides,1);   cPeps = repmat({'',0},numPeps,1);
pepsIdx = zeros(numPeps,1);

idxNE = find(strcmp(peptides,'')==0);   if isempty(idxNE) return; end;

pepsNE = peptides(idxNE,:);   szNE = length(idxNE);
[pepsS,idxS]=sort(pepsNE);   idxNE = idxNE(idxS);
idxNewPep = [1+find(strcmp(pepsS(1:szNE-1),pepsS(2:szNE))==0)];
rangeSame = [[1;idxNewPep] [idxNewPep-1;szNE]];   szRange = size(rangeSame,1);

pepsU = peptides(idxNE(rangeSame(:,1)));
for i=1:szRange
    idx=idxNE(rangeSame(i,1):rangeSame(i,2));   curPep=pepsU{i};   pepLen=length(curPep);   result={'',0};
    
    % Build the unmodified version of a peptide - remove [123] and *,#,others. Count # mods.
    if pepLen>0
        idxS = find(curPep=='[');   idxE = find(curPep==']');
        if size(idxS,2)==size(idxE,2)
            numMods = size(idxS,2);  toKeep = ones(1,pepLen);
            for m=1:numMods toKeep(idxS(m):idxE(m))=0; end;
            
            newPep = curPep(find(toKeep==1));   newLen = length(newPep);
            posOk = find(newPep>='A' | newPep<='Y');   szOk = length(posOk);
            if szOk~=newLen
                newPep = newPep(posOk);   numMods=numMods+newLen-szOk;
            end;
            result{1} = newPep;   result{2} = numMods;
        else fprintf(1,'Peptide %d has an incorrect format: %s\n',i,curPep); end;
    end;
    
    cPeps(idx,:) = repmat(result,rangeSame(i,2)-rangeSame(i,1)+1,1);
end
