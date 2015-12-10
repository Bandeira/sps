function clustersIdx = save_clusters(filename, seqSpecs, vSets, eSets, aligns, tolerance)
% function clustersIdx = save_clusters(filename, seqSpecs, vSets, eSets, aligns, tolerance)
%
%  Save clusters information for c++ code
%
%  seqSpecs - as returned by denovoMA, 3 col cell: consensus PRMs (col 1), endpoints (col 2), parent mass (col 3)
%
%  clustersIdx - lists the cluster/subcluster per saved entry in the clusters files (filename)
%

fid=fopen(filename,'w');  if fid<=0 fprintf(1,'ERROR opening %s!\n',filename); return; end;

numSubClusters = 0; for i=1:size(vSets,1) numSubClusters=numSubClusters+size(vSets{i,1},1); end;
clustersIdx = zeros(numSubClusters,2);   curCluster = 1;
fprintf(fid,'%.0d\n',numSubClusters);
for i=1:size(vSets,1)
    for j=1:size(vSets{i,1},1)   % Treat each subcluster as a separate cluster
        clustersIdx(curCluster,:) = [i j]; curCluster = curCluster+1;
        
        % Save left/right shifts per spectrum in cluster
        fprintf(fid,'%.0d\n',size(vSets{i,1}{j},2));
        shifts1 = getAbsoluteShifts(aligns(eSets{i}{j,1},:),tolerance,0);
        shifts2 = getAbsoluteShifts(aligns(eSets{i}{j,2},:),tolerance,0);
        for k=1:size(vSets{i,1}{j},2)
            fprintf(fid,'%.0d %.2f %.2f\n',vSets{i,1}{j}(k), shifts1(vSets{i,1}{j}(k)), shifts2(vSets{i,1}{j}(k)));
        end;
        
        % Save consensus spectrum information (if available)
        if ~isempty(seqSpecs)
            % Consensus PRMs
            fprintf(fid,'%.0d %.1f\n',size(seqSpecs{i}{j}{1},1),seqSpecs{i}{j}{3});  % Num. peaks ; Parent mass
            for k=1:size(seqSpecs{i}{j}{1},1) fprintf(fid,'%.1f %.1f\n',seqSpecs{i}{j}{1}(k,1:2)); end;
            
            % Endpoints
            fprintf(fid,'%.0d\n',size(seqSpecs{i}{j}{2},1));   % Num. endpoints
            for k=1:size(seqSpecs{i}{j}{2},1) fprintf(fid,'%.1f %.1f\n',seqSpecs{i}{j}{2}(k,1:2)); end;
        else
            fprintf(fid,'0 0\n0\n');
        end
    end
end

fclose(fid);
