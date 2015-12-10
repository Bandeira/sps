function [vSets, eSets] = assembleE(aligns)
% function [vSets, eSets] = assembleE(aligns)
%
% Find all the connected components in an alignment graph. 
%
%   vSets - sets of vertices in each connected component
%   eSets - sets of edges in each connected component (as indices of aligns)

% numPeptides = size(contigs,1);
numPeptides = size(unique(aligns(:,[1 2])),1);

vSets = cell(numPeptides,2);

contigsIdx = zeros(max(max(aligns(:,[1 2]))),1);
numContigSets = 0;

%errCount=0;

szAligns = size(aligns,1);
toProcess = ones(szAligns,1); idxAligns = [1:szAligns]';   tensLeft = floor(szAligns/10000);
curTime = 0; tic;
% while size(idxAligns,1)>0
%     [score, idx] = max(aligns(idxAligns,4));   idx = idxAligns(idx(1));
for idx=1:szAligns
    
    i = aligns(idx,1);
    j = aligns(idx,2);
%     shift = aligns(idx,3);
    
    if contigsIdx(i)>0 & contigsIdx(j)>0 & contigsIdx(i)==contigsIdx(j) % ignore edges between vertices already in the same cluster - MAKES eSets INCOMPLETE !!!!!!!!
        continue; 
    end;  
    
%     relatedEdges = find((aligns(:,1)==j & aligns(:,2)==i) | (aligns(:,1)==i & aligns(:,2)==j))';

%     idxI = findContig(i,vSets);      % 04/9/14 - Changed for efficiency purposes
%     idxJ = findContig(j,vSets);
    idxI = contigsIdx(i);
    idxJ = contigsIdx(j);
    
    if idxI==0 & idxJ==0 numContigSets=numContigSets+1; vSets{numContigSets,1} = [i j];   contigsIdx([i j])=numContigSets;   
    else if idxI==0 vSets{idxJ,1} = [vSets{idxJ,1} i];  contigsIdx(i)=idxJ;  
        else if idxJ==0 vSets{idxI,1} = [vSets{idxI,1} j];  contigsIdx(j)=idxI;
            else if idxI~=idxJ 
                    vSets{idxI,1} = unique([vSets{idxI,1} vSets{idxJ,1}]); 
                    contigsIdx(vSets{idxJ,1}) = idxI;
                    if idxJ<numContigSets
                        vSets{idxJ,1} = vSets{numContigSets,1}; 
                        contigsIdx(vSets{numContigSets,1}) = idxJ;
                    end
                    numContigSets=numContigSets-1;
                end;
            end;
        end;
    end;

    if max(contigsIdx)>numContigSets 
        a=1;
    end
    
%    fprintf(1,'New contig from %d to %d (%d,%d), score = %.2f confidence = %2.3f, aligns is now size %d\n',i,j,posI,posJ,aligns(idx(1),4),aligns(idx(1),5),size(aligns,1));
    

%     toProcess(relatedEdges)=0;
%     idxAligns = find(toProcess==1);
% %     idxAligns = setdiff(idxAligns, relatedEdges)';   % 'discard' all edges whose end vertices were just processed
   
    curTensLeft=floor(size(idxAligns,1)/10000); 
    if curTensLeft<tensLeft 
        t=toc; curTime=curTime+t; tic;
        fprintf(1,'%d aligns left, current 10000 took %.1f secs, ETA: %.1f secs...\n',size(idxAligns,1),t,curTensLeft*curTime/10000); tensLeft=curTensLeft; 
    end;
end; % while
vSets = vSets(1:numContigSets,:);
edgesSet = contigsIdx(aligns(:,1));   eSets = cell(size(vSets,1),1);
for i=1:size(vSets,1)  eSets{i} = find(edgesSet==i)'; end;

function idx = findContig(vertex, vSets)

idx=0;
for c=1:size(vSets,1)
    if ismember(vertex, vSets{c,1}) idx=c; return; end;
end;
