function [pairs, aligns, ratios, means] = merge_batch_asp2(resultFiles, ratioFiles, meanFiles)
%
%  Merges results files from batch_asp. All 3 parameter arguments are cell vectors with 1 col.
%
%  Similar to merge_batch_asp.m but eliminates eventually repeated pairs by keeping only the
%    alignment with the highest aggregated score (sum of the match scores for both spectra).
%

numFiles = size(resultFiles,1);
pairs = []; aligns = []; ratios = []; meansTmp = cell(1,numFiles);   meansCorrections = [];
for i=1:numFiles
    [pairsCur, alignsCur] = load_resultsASP(resultFiles{i},';');
    
    % Ensure that pairs has the format [lower index, higher index] - Added Sep.8,2005
    idx = find(pairsCur(:,1)>pairsCur(:,2));   
    pairsCur(idx,:) = pairsCur(idx,[2 1]);    alignsCur(idx,:) = [single(-double(alignsCur(idx,1))) alignsCur(idx,[3 2])];
    
    [common, idxOld, idxCur] = intersect(pairs,pairsCur,'rows');
    if isempty(common)
        pairs = [pairs; pairsCur];   aligns = [aligns; alignsCur];
        ratios = [ratios; load(ratioFiles{i})];
        meansTmp{i} = load(meanFiles{i});
    else
        ratiosCur = load(ratioFiles{i});   meansTmp{i} = load(meanFiles{i});
        idxNew = setdiff([1:size(pairsCur,1)]',idxCur);
        pairs = [pairs; pairsCur(idxNew,:)];   aligns = [aligns; alignsCur(idxNew,:)];   ratios = [ratios; ratiosCur(idxNew,:)];
        
        % Select between repeated pairs and correct means
        for j=1:size(common,1)
            if (double(aligns(idxOld(j),2)) + double(aligns(idxOld(j),3))) >= (double(alignsCur(idxCur(j),2)) + double(alignsCur(idxCur(j),3)))
                scoresToRemove = double(alignsCur(idxCur(j),2:3));
            else
                scoresToRemove = double(aligns(idxOld(j),2:3));
                aligns(idxOld(j),:) = alignsCur(idxCur(j),:);
                ratios(idxOld(j),:) = ratiosCur(idxCur(j),:);
            end

            spec1 = pairs(idxOld(j),1);   spec2 = pairs(idxOld(j),2);
            meansCorrections(spec1,:) = meansCorrections(spec1,:) + [scoresToRemove(1,1) 1];
            meansCorrections(spec2,:) = meansCorrections(spec2,:) + [scoresToRemove(1,2) 1];
        end
    end
    
    if isempty(meansCorrections) meansCorrections=zeros(size(meansTmp{1})); end;
    fprintf(1,'Finished processing %s\n',resultFiles{i});
end;

% Merge all the means
proportions = zeros(size(meansTmp{1},1),numFiles);
for i=1:numFiles  proportions(:,i) = meansTmp{i}(:,2);  end;  % Get the counts
counts = sum(proportions')';
for i=1:numFiles  proportions(:,i) = proportions(:,i) ./ counts;  end;

means = zeros(size(proportions,1),1);
for i=1:numFiles  means = means + proportions(:,i).*meansTmp{i}(:,1);  end;

% Correct means for the repeated pairs
proportions = counts./(counts-meansCorrections(:,2));   values = meansCorrections(:,1)./counts;
means = proportions.*(means-values);
