function [vSets_d, eSets_d] = addVerticesSP(vSets_allV, eSets_allV, vSets_d, eSets_d, aligns, shiftTolerance)
%
%  Same Peptide vertices are not adequately merged by filterDirection3. This function considers the
%    vertices left out by filterDirection3 and, if they're the same as some spectrum in vSets_d (both shifts are zero)
%    then add the vertex to vSets_d
%
%  vSets_allV, eSets_allV - as returned by assembleE
%  vSets_d, eSets_d - as returned by filterDirection3
%  shiftTolerance   - Tolerance to consider a shift as a zero shift (e.g. 2 Da)
%

fprintf(1,'Warning: Keeping edge indices as subindices of eSets_all instead of direct indices into aligns\n');

numSets = size(vSets_allV,1);
for i=1:numSets
    if isempty(vSets_d{i}) fprintf(1,'Warning: Cluster %d has no subclusters. Skipping...\n',i); continue; end;
    if size(vSets_allV{i},2)==size([vSets_d{i}{:}],2) continue; end;
    
    allV = vSets_allV{i};         numVertices = size(allV,2);
    alignsV = aligns(eSets_allV{i},:);
    
    renameVec = zeros(max(allV),1);   % Rename the vertices to 1:numVertices to save memory (on the adjacency matrix)
    renameVec(allV) = [1:size(allV,2)];
    alignsV(:,1:2) = sort(renameVec(alignsV(:,1:2))')';
    
    adjacencySP = zeros(numVertices, numVertices);  % [i,j]==1 if spectra i and j have only zero shifts between them
    uniqueEdges = unique(alignsV(:,1:2),'rows');
    for j=1:size(uniqueEdges,1)
        shifts = abs(alignsV(find(alignsV(:,1)==uniqueEdges(j,1) & alignsV(:,2)==uniqueEdges(j,2)),3));
        if isempty(find(shifts>shiftTolerance)) 
            adjacencySP(uniqueEdges(j,1),uniqueEdges(j,2))=1; adjacencySP(uniqueEdges(j,2),uniqueEdges(j,1))=1; 
        end;
    end
    
    vertexMembership = zeros(numVertices,1);        % For each vertex, indicates the subcluster it belongs to
    szSubClusters    = zeros(size(vSets_d{i},1),1); % Get the size of each subcluster
    for j=1:size(vSets_d{i},1)
        vertexMembership(renameVec(vSets_d{i}{j})) = j;    
        szSubClusters(j)=size(vSets_d{i}{j},2);
    end
    
    newV = renameVec(setdiff(allV',[vSets_d{i}{:}]'));
    lastCount=0;
    while ~isempty(newV) & size(newV,1)~=lastCount                   % Repeat the whole process until every connected spectrum is placed somewhere
        lastCount = size(newV,1);   processed = zeros(size(newV,1),1);
        for j=1:size(newV,1)
            adj = find(adjacencySP(newV(j),:)==1)';   if isempty(adj) continue; end;
            memb = vertexMembership(adj);             if max(memb)==0 continue; end;

            memb = memb(find(memb>0));
            if max(size(unique(memb)))>1 
                counts = hist(memb,[0:max(memb)]);
                memb = memb( find( counts(memb+1) == max(counts) ) );   % Consider only the clusters to which newV(i) is maximally connected to
            end
            
            if max(size(unique(memb)))>1 
                [foo, idxS] = sort(szSubClusters(memb));  % Find largest cluster
                memb = memb(idxS(size(idxS,1)));
            end
            
            memb = memb(1);  % Just to make sure that memb has only a single value
            
            % Merge new vertex with selected subcluster memb
            vSets_d{i}{memb} = [vSets_d{i}{memb} allV(newV(j))];
            idx = find(vertexMembership(adj)==memb);      % Find other SP vertices in the same subcluster and collect edges to them
            for k=1:size(idx,1)
                edgeIdx = find(alignsV(:,1)==min(newV(j),adj(idx(k))) & alignsV(:,2)==max(newV(j),adj(idx(k))));
                if ~isempty(edgeIdx)
                    if size(edgeIdx,1)<2 fprintf(1,'ERROR: Cluster %d - there is only one shift between vertices %d and %d !?!\n',i,allV([newV(j) adj(idx(k))])); 
                    else
%                         eSets_d{i}{memb,1}=[eSets_d{i}{memb,1} eSets_allV{i}(edgeIdx(1))];
%                         eSets_d{i}{memb,2}=[eSets_d{i}{memb,2} eSets_allV{i}(edgeIdx(2))];
% Keep edge indices as subindices of eSets_all instead of direct indices into aligns
                        eSets_d{i}{memb,1}=[eSets_d{i}{memb,1} edgeIdx(1)];
                        eSets_d{i}{memb,2}=[eSets_d{i}{memb,2} edgeIdx(2)];
                    end
                end
            end
            
            processed(j)=1;
        end
        newV = newV(find(processed==0));
    end
end