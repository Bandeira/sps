function [vSets, eSets,trigVerts] = filterDirection3(contigs, aligns, edgeScores, tolerance)
%
%  filterDirection3: Same generic purpose as versions 1 and 2 - find subsets of coherent edges. Differences are:
%     - Intended for usage on large sets of edges, including sets with incorrect pairwise alignments
%     - Symmultaneously constructs the symmetric coherent edgesets and actively avoids merging them (by avoiding edges that would merge the sets)
%     - Selects edges to use by decreasing order of edgeScores (highest scoring edge first)
%     - Symmultaneous construction of different connected components
%
%  vSets, eSets - lists of vertices and edges for the connected components
%  trigVerts - Lists of vertices of the leftover unprocessed triangles (not enough information to determine merge direction)
% 

% NOTE: Correct bug in getMatchRatios that inverts the ratios assuming that the symmetric of v1->v2 is always v2->v1

% List of all the important variables in this function and their purpose
%
% vertices   - used to rename the vertices (to save memory in the internal processing variables)
% edges      - adjacency matrix; each entry [i,j] contains the indices of the _two_ edges between v_i and v_j (no zero shift redundancy)
% adjacency  - adjacency list; useful when intersecting adjacencies
%
% trigSets   - line i contains the sets of edges that define the two symmetric versions of the i-th triangle (one in col 1 and the other in col 2)
% trigVerts  - line i contains the sorted set of vertices that define the i-th valid triangle 
% trigScores - entry i contains the score of the i-th triangle
%
% vSets      - vertices of the connected components
% eSets      - edges of the connected components; column 1 contains edgesets in one direction and column 2 contains the symmetric direction
% membership - indicates which component each vertex belongs to
% edgeMemberships - indicates which eSets direction an edge belongs to (either 1 or 2, 0 for non-assigned edges)
% 

DEBUG=1;

szAligns = size(aligns,1);
if size(aligns,2)~=5 fprintf(1,'Error: filterDirection3 assumes that the parameter aligns has 5 columns. Exiting...\n'); return; end;

% Rename vertices to save memory - REMEMBER TO RENAME THEM BACK AT THE END !!!!!!!!!
vertices = unique([aligns(:,1); aligns(:,2)]);    numVertices = size(vertices,1);
renameVec = zeros(max(vertices),1);               renameVec(vertices) = [1:numVertices];
aligns = [renameVec(aligns(:,1:2)) aligns(:,3:size(aligns,2))];
contigs = contigs(vertices,:);


% Reorder all shifts to [min vertex idx, max vertex idx, positive/negative shift]
for i=1:szAligns
    if aligns(i,1)>aligns(i,2) 
        aligns(i,1:2) = [aligns(i,2) aligns(i,1)];
        if aligns(i,3)>0 aligns(i,3) = -aligns(i,3); end;
    end
end
[foo, idxSaligns] = sortrows(aligns,1:3);   aligns=aligns(idxSaligns,:);   edgeScores=edgeScores(idxSaligns,:);   clear foo;


% Remove double zero shift edges + Mark edges with equal symmetric shifts (both equal to middle)
edges = cell(numVertices,numVertices);  % edges{v1,v2} - contains edges from v1 to v2 (positive or negative)
adjacency  = cell(numVertices,1);       % lists of adjacent vertices (indexed per vertex)
alignmentScores = zeros(size(unique(aligns(:,1:2),'rows'),1),3);   % Contains the scores per pairwise alignment (not just per edge)
alignmentScoresIdx = 1;
i=1;
while i<szAligns
    v1 = aligns(i,1);   v2 = aligns(i,2);
    if aligns(i+1,1)~=v1 | aligns(i+1,2)~=v2 i=i+1; continue; end;  % Symmetric shift is missing
    alignmentScores(alignmentScoresIdx,1:2)=[v1 v2];
    if i+2<=szAligns & aligns(i+2,1)==v1 & aligns(i+2,2)==v2
        % This can only happen when one of the shifts is non-zero and the other is zero (therefore appearing twice)
        alignmentScores(alignmentScoresIdx,:)=max(edgeScores(i:i+2));   
        % The edge in the middle (i+1) should always be a zero shift - only need to find the non-zero shift (either at i or i+2)
        if abs(aligns(i,3))<=tolerance edges{v1,v2}=[i+1 i+2]; else edges{v1,v2}=[i i+1]; end;
        i=i+3;
    else 
        edges{v1,v2} = [i i+1];  alignmentScores(alignmentScoresIdx,:)=max(edgeScores(i:i+1));   i=i+2;
    end
    
    adjacency{v1} = [adjacency{v1}; v2];
    adjacency{v2} = [adjacency{v2}; v1];
end

%  - What are the edges with equal symmetric shifts?  Used when adding extra triangles defined by vertices already in the component
idx = [edges{:,:}];
[idxGood, idxZero, idxEqualLR] = separateEdges( aligns(idx,:), tolerance );
degenerateEdges = zeros(szAligns,1);   degenerateEdges(idx([idxZero;idxEqualLR])) = 1;

%  - What are the possible triangles and their scores?
trigSets = {};                    % Keeps track of symmetric triangles; each position contains a coherent triangle
trigVerts = [];                   % Array with sorted list of the triangles' vertices (for fast location of triangle edges in trigSets)
trigScores = [];                  % Vector of triangle scores (sum of edges' scores)
for v1=1:numVertices
    adj_v1 = adjacency{v1};   adj_v1 = adj_v1(find(adj_v1>v1));   % Triangles are considered in increasing order of vertex indices only
    for j=1:size(adj_v1,1)
        v2 = adj_v1(j);
        adj_v1v2 = intersect(adjacency{v1},adjacency{v2});
        adj_v1v2 = adj_v1v2(find(adj_v1v2>v2));
        for k=1:size(adj_v1v2,1)
            v3 = adj_v1v2(k);
            e12 = [aligns(edges{v1,v2},:) edges{v1,v2}'];    % Keep track of the edges' indices in the last column (col 6)
            e13 = [aligns(edges{v1,v3},:) edges{v1,v3}'];
            e23 = [aligns(edges{v2,v3},:) edges{v2,v3}'];
            
            [t1,t2] = getTriangles(e12,e13,e23,4*tolerance);
            if isempty(t1) continue; end;
            trigSets = [trigSets; {t1} {t2}];    trigVerts = [trigVerts; [v1 v2 v3]];   trigScores = [trigScores; sum(edgeScores(t1(:,6)))];
        end
    end
end
numTrigs = size(trigScores,1);
[foo,idxS] = sort(trigScores);   idxS=idxS(numTrigs:-1:1);   trigScores = trigScores(idxS);
trigSets=trigSets(idxS,:);       trigVerts = trigVerts(idxS,:);
fprintf(1,'Got %d triangles\n',numTrigs);

%
%  ---------------------------------------------------
%   Phase 2 - Build the coherent connected components
%  ---------------------------------------------------
%

vSets = {};    eSets = {};    numComps=0;
membership = zeros(numVertices,1);    % Keep track of what component each vertex belongs to
edgeMembership = zeros(szAligns,1);   % Keep track of which side (dir/sym) of a component an edge belongs to

% numIterations = 0;   tic;
while numTrigs>0
    i=1;
    while i<=numTrigs
        comps = membership(trigVerts(i,:));    [compsS, idxS] = sort(comps);
        if compsS(3)==0    % New component
            vSets = [vSets; {[trigVerts(i,:)]}];
            eSets = [eSets; {trigSets{i,1}(:,6)'} {trigSets{i,2}(:,6)'}];
            numComps = numComps+1;
            membership(trigVerts(i,:))=numComps;
            edgeMembership(eSets{numComps,1})=1;
            edgeMembership(eSets{numComps,2})=2;
            processedTrigs = [i];   processedComp = numComps;
            break;
        end
        
        if compsS(1)==0 & compsS(2)==compsS(3)    % Adding a single vertex to a component
            c = compsS(3);   newV = trigVerts(i,idxS(1));

            eMemb = edgeMembership(trigSets{i,1}(:,6));  setDir1 = max(eMemb);   newEdges1 = find(eMemb==0);   % setDir should only be 1 or 2
            if degenerateEdges(trigSets{i,1}(find(eMemb>0),6))==1 i=i+1; continue; end;   % Can't add new vertex if the edge in the component is degenerate
            
            eMemb = edgeMembership(trigSets{i,2}(:,6));  setDir2 = max(eMemb);   newEdges2 = find(eMemb==0);   % and the two setDirs should always be different
            if setDir1==setDir2 | size(newEdges1,1)~=2 | size(newEdges2,1)~=2
                fprintf(1,'ERROR: directions should be unique!!\n');
                vSets = {};   eSets = {};   return;
%                 a=1/0;
            end
            eSets{c,setDir1}=[eSets{c,setDir1} trigSets{i,1}(newEdges1,6)'];   edgeMembership(trigSets{i,1}(newEdges1,6))=setDir1;
            eSets{c,setDir2}=[eSets{c,setDir2} trigSets{i,2}(newEdges2,6)'];   edgeMembership(trigSets{i,2}(newEdges2,6))=setDir2;

            vSets{c} = [vSets{c} newV];
            membership(newV) = c;

            processedTrigs = [i];   processedComp = c;
        end
        
        if compsS(1)>0 & (compsS(1)==compsS(2) | compsS(2)==compsS(3))   % Attempting merge of two components
            c2 = compsS(2);   v2=trigVerts(i,idxS(2));   % position 2 must contain the index of the component with 2 selected vertices
            if c2==compsS(1)
                c1 = compsS(3);  v1=trigVerts(i,idxS(1));  vNew=trigVerts(i,idxS(3));
            else
                c1 = compsS(1);  v1=trigVerts(i,idxS(3));  vNew=trigVerts(i,idxS(1));
            end
            if degenerateEdges(find(aligns(:,1)==min(v1,v2) & aligns(:,2)==max(v1,v2)))==1 i=i+1; continue; end; % Degenerate edges won't determine direction
            
            setBoth = intersect(vSets{c1}',adjacency{vNew}); % Set of vertices eligible to be the second vertex in c1
            if isempty(setBoth) i=i+1; continue; end;        % Try next best triangle
            
            hyp1 = intersect(adjacency{v1},setBoth);   hyp2 = intersect(adjacency{v2},setBoth);
            hyp = [v1*ones(size(hyp1)) hyp1; v2*ones(size(hyp2)) hyp2];
            for j=1:size(hyp,1)
                if degenerateEdges(find(aligns(:,1)==min(vNew,hyp(j,2)) & aligns(:,2)==max(vNew,hyp(j,2))))==1 continue; end; % Degenerate edges won't determine direction
                trigTry = sort([hyp(j,1) vNew hyp(j,2)]);
                trigTryIdx = find(trigVerts(:,1)==trigTry(1) & trigVerts(:,2)==trigTry(2) & trigVerts(:,3)==trigTry(3));
                if isempty(trigTryIdx) continue; end;

                % Are the trigs coherent?
                [trigMergeDir,trigMergeSym] = intersectEdgesets(trigSets{i,1},trigSets{i,2},trigSets{trigTryIdx,1},trigSets{trigTryIdx,2},tolerance);
                if isempty(trigMergeDir) continue; end;
                
                % What set of edges matches the triangles for both components?
                idx2 = find(trigMergeDir(:,1)==min(v1,v2) & trigMergeDir(:,2)==max(v1,v2));  % Get the membership of the reference edges
                idx1 = find(trigMergeDir(:,1)==min(vNew,hyp(j,2)) & trigMergeDir(:,2)==max(vNew,hyp(j,2)));
                if c2<c1  
                    mergeDirs = edgeMembership(trigMergeDir([idx2 idx1],6))';   % And determine direction of merge based on those reference edges
                    processedComp = c2; otherComp = c1;
                else
                    mergeDirs = edgeMembership(trigMergeDir([idx1 idx2],6))';
                    processedComp = c1; otherComp = c2;
                end
                mergeDirs = [mergeDirs; mod(mergeDirs,2)+1];   % Note: second line is just complement of first line 1,2<->2,1 and 1,1<->2,2
                
                newEdgesD = find(edgeMembership(trigMergeDir(:,6))==0);
                newEdgesS = find(edgeMembership(trigMergeSym(:,6))==0);
                
                % Merge the components
                vSets{processedComp} = [vSets{processedComp} vSets{otherComp}];
                eSets{processedComp,mergeDirs(1,1)} = [eSets{processedComp,mergeDirs(1,1)} eSets{otherComp,mergeDirs(1,2)} trigMergeDir(newEdgesD,6)'];
                eSets{processedComp,mergeDirs(2,1)} = [eSets{processedComp,mergeDirs(2,1)} eSets{otherComp,mergeDirs(2,2)} trigMergeSym(newEdgesS,6)'];
                membership(vSets{otherComp}) = processedComp;                       % Change vertex memberships according to new merged component
                edgeMembership(eSets{otherComp,mergeDirs(1,2)}) = mergeDirs(1,1);   % and edgeMemberships too
                edgeMembership(trigMergeDir(newEdgesD,6)) = mergeDirs(1,1);
                edgeMembership(eSets{otherComp,mergeDirs(2,2)}) = mergeDirs(2,1);
                edgeMembership(trigMergeSym(newEdgesS,6)) = mergeDirs(2,1);
                
                if otherComp<numComps  
                    vSets{otherComp} = vSets{numComps};   eSets(otherComp,:) = eSets(numComps,:);   membership(vSets{numComps}) = otherComp;
                end
                numComps=numComps-1;   vSets=vSets(1:numComps,:);   eSets = eSets(1:numComps,:);
                
                processedTrigs = [i trigTryIdx];
                break;
            end
        end
        
        if DEBUG==1
            idxNonDeg = eSets{processedComp,1}(find(degenerateEdges(eSets{processedComp,1})==0));
            if size(unique(edgeMembership(eSets{processedComp,1})),1)==2
                fprintf(1,'Found component with incoherent directions!\n');
                a=1/0;
            end
            idxNonDeg = eSets{processedComp,2}(find(degenerateEdges(eSets{processedComp,2})==0));
            if size(unique(edgeMembership(eSets{processedComp,2})),1)==2
                fprintf(1,'Found component with incoherent directions!\n');
                a=1/0;
            end
        end
        
        % Add other elegible triangles to the merged components
        if ~isempty(processedTrigs)
            % Get all triangles from processedComp still in trigVerts
            trigMship = membership(trigVerts);   if size(trigVerts,1)==1  trigMship=trigMship'; end;
            extraTrigs = setdiff(find(trigMship(:,1)==processedComp & trigMship(:,2)==processedComp & trigMship(:,3)==processedComp),processedTrigs');
            
            for j=1:size(extraTrigs,1)
                degOk = find(degenerateEdges(trigSets{extraTrigs(j),1}(:,6))==0);  % Only non-degenerate alignment edges can be used to determine direction of merge
                eMemb1 = edgeMembership(trigSets{extraTrigs(j),1}(:,6));  setDir1 = max(eMemb1(degOk));   newEdges1 = find(eMemb1==0);   % setDir should only be 1 or 2
                degOk = find(degenerateEdges(trigSets{extraTrigs(j),2}(:,6))==0);
                eMemb2 = edgeMembership(trigSets{extraTrigs(j),2}(:,6));  setDir2 = max(eMemb2(degOk));   newEdges2 = find(eMemb2==0);   % and the two setDirs should always be different

                if (isempty(newEdges1) | isempty(newEdges2)) | (isempty(setDir1) | isempty(setDir2)) | (setDir1==0 | setDir2==0)
                % first () : In this case there are no new edges so skip remaining steps
                % second (): there are no non-degenerate edges in the triangle
                % third () : ==0 means that only the non-assigned edges are non-degenerate
                   continue; 
               end;   

               if setDir1==setDir2   
                   fprintf(1,'ERROR: directions should be unique when adding extra triangle edges to the processed component!! Trig = (%d,%d,%d)\n',vertices(trigVerts(extraTrigs(j),:)));
%                    fprintf(1,'Terminating with %d triangles left unprocessed.\n',numTrigs-size(processedTrigs,2)-j+1);
                    i=numTrigs+1; % FORCE EXIT AT END OF LOOP
                    break;
%                    vSets = {};   eSets = {};   return;
               end

               eSets{processedComp,setDir1}=[eSets{processedComp,setDir1} trigSets{extraTrigs(j),1}(newEdges1,6)'];   edgeMembership(trigSets{extraTrigs(j),1}(newEdges1,6))=setDir1;
               eSets{processedComp,setDir2}=[eSets{processedComp,setDir2} trigSets{extraTrigs(j),2}(newEdges2,6)'];   edgeMembership(trigSets{extraTrigs(j),2}(newEdges2,6))=setDir2;
            end
            
            processedTrigs = [processedTrigs extraTrigs'];
            % Note: by construction, merging never creates a cycle because the sets A & B are disjoint prior to merging
            %         Only this addition of extra triangles may create cycles in the graph (this IS the desired behaviour)
            break;
        end
        
        i=i+1;
    end
    if i>numTrigs
        fprintf(1,'Terminating with %d triangles left unprocessed\n',numTrigs);
        break;
    end
    
    if DEBUG==1
        idxNonDeg = eSets{processedComp,1}(find(degenerateEdges(eSets{processedComp,1})==0));
        if size(unique(edgeMembership(eSets{processedComp,1})),1)==2
            fprintf(1,'Found component with incoherent directions!\n');
            a=1/0;
        end
        idxNonDeg = eSets{processedComp,2}(find(degenerateEdges(eSets{processedComp,2})==0));
        if size(unique(edgeMembership(eSets{processedComp,2})),1)==2
            fprintf(1,'Found component with incoherent directions!\n');
            a=1/0;
        end
    end
    
    kept = setdiff([1:numTrigs]',processedTrigs');
    trigSets = trigSets(kept,:);    trigVerts = trigVerts(kept,:);   trigScores = trigScores(kept);
    numTrigs = size(trigSets,1);    processedTrigs = [];
    
%     numIterations = numIterations + 1;
%     fprintf(1,'Finished iteration %d in %.1f seconds. Triangles left: %d\n',numIterations, toc, numTrigs); pause(.01); tic;
end

% vertices were renamed to save memory; just renaming them back at the end
for i=1:size(vSets,1)  vSets{i} = vertices(vSets{i})'; end;
% aligns was sorted so convert edge indices to pre-sorting indices
for i=1:size(eSets,1)  eSets{i,1} = idxSaligns(eSets{i,1})';   eSets{i,2} = idxSaligns(eSets{i,2})'; end;



function [setRdir,setRsym] = intersectEdgesets(set1dir,set1sym,set2dir,set2sym,tolerance)
%
%  Coherently intersect edgesets - some edge must uniquely belong to a combination of set1dir/sym and set2dir/sym
%
%  NOTE: There can be only one common edge when intersecting two triangles and also when intersecting
%          two merged triangles with two disjoint coherent sets

% fprintf(1,'WARNING: Degenerate shifts cannot be used to determine correct direction of alignment\n');

match = repmat(set1dir([1 1 1 2 2 2 3 3 3]',1:3),2,1) - [repmat(set2dir(:,1:3),3,1); repmat(set2sym(:,1:3),3,1)];
idx = find(match(:,1)==0 & match(:,2)==0 & abs(match(:,3))<=tolerance);

if isempty(idx) | max(size(idx))>1
    % NOTE: Any two different triangles can have at most one edge in common. If this single edge from set1 matches
    %        both an edge from set2dir and an edge from set2sym then the edge from set1 is degenerate and cannot
    %        be used to coherently merge the triangles.
    setRdir=[]; setRsym=[]; return;
end

if idx>9   % Switch set2sym to match set1dir
    idx = idx-9;
    tmp = set2dir;   set2dir=set2sym;   set2sym=tmp;
end
idx = mod(idx,3);  if idx==0 idx=3; end;

setRdir=[set1dir; set2dir(1:idx-1,:); set2dir(idx+1:3,:)];
idx2=find(set2sym(:,1)==set2dir(idx,1) & set2sym(:,2)==set2dir(idx,2)); setRsym=[set1sym; set2sym(1:idx2-1,:); set2sym(idx2+1:3,:)];
