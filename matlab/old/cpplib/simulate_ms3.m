function [idxBase, sets] = simulate_ms3(pos, maxSubSpecs, clusterFN)
% function [idxBase, sets] = simulate_ms3(pos, maxSubSpecs, clusterFN)
%
%  Simulates MS3 datasets from prefix/suffix spectral pairs. Each peptide that 
%   contains at least one other prefix/suffix peptide (according to pos) will generate one
%   MS3-set.
%
%  pos(i,:) - start/end peptide positions on protein (cols 1/2)
%  maxSubSpecs - If not zero then each MS3-set will contain at most maxSubSpecs child spectra
%                 which are selected randomly from all the possible child spectra
%  clusterFN - name of the file where the sets should be saved to
%
%  idxBase   - index of the 'parent MS3' spectra output to clusterFN
%  sets      - the simulated MS3 sets. Each cell entry is [parent_spectrum <child spectra>]
% 

% [foo,idxS] = sortrows(pos);   pos = pos(idxS,:);   
numSpecs = size(pos,1);

sets = cell(numSpecs,1);   idxBase = zeros(numSpecs,1);   setIdx=1;
for specIdx=1:numSpecs
    idxChildren = find((pos(:,1)==pos(specIdx,1) & pos(:,2)<pos(specIdx,2)) | (pos(:,1)>pos(specIdx,1) & pos(:,2)==pos(specIdx,2)) );
    if isempty(idxChildren) continue; end;

    if size(idxChildren,1)>maxSubSpecs idxChildren=idxChildren(randperm(size(idxChildren,1))); idxChildren=idxChildren(1:maxSubSpecs); end;
    sets{setIdx} = [specIdx idxChildren'];   idxBase(setIdx) = specIdx;   setIdx=setIdx+1;
end;

sets = sets(1:setIdx-1,:);   idxBase = idxBase(1:setIdx-1,:);   numSets = size(sets,1);   
fprintf(1,'Generated %d MS3 sets.',numSets); if numSets==0 return; end;
fid = fopen(clusterFN,'w'); if fid<=0 fprintf(1,'ERROR opening %s!\n',clusterFN); return; end;
fprintf(fid,'%.0f\n',numSets);
for setIdx=1:numSets
    fprintf(fid,'%.0f\n',size(sets{setIdx},2));
    for specIdx=1:size(sets{setIdx},2) fprintf(fid,'%.0f 0 0\n',sets{setIdx}(specIdx)-1); end;  % Zero-based spectrum indices
    fprintf(fid,'0 0\n');  % No PRM entries
    fprintf(fid,'0\n');    % No endpoint entries
end
fclose(fid);

