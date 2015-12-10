function [pairs, aligns, ratios, means, stddevs] = merge_batch_pa3(resultFiles, ratioFiles, meanFiles, varsFiles)
% function [pairs, aligns, ratios, means, stddevs] = merge_batch_pa3(resultFiles, ratioFiles, meanFiles, varsFiles)
%
%  Merges results files from batch_asp. All 3 parameter arguments are cell vectors with 1 col.
%
%  Similar to merge_batch_asp3.m but for pa aligns.
%

numFiles = size(resultFiles,1);
pairs = []; aligns = []; ratios = []; meansTmp = cell(1,numFiles);   meansCorrections = [];   varsTmp = cell(1,numFiles);   varsCorrections = [];
for i=1:numFiles
    [pairsCur, alignsCur] = load_resultsPA(resultFiles{i},';');
    meansTmp{i} = load(meanFiles{i});   if i==1 meansCorrections=zeros(size(meansTmp{1},1),2); end;
    varsTmp{i} = load(varsFiles{i});    if i==1 varsCorrections=zeros(size(meansTmp{1},1),1); end;
    
    % Ensure that pairs has the format [lower index, higher index] - Added Sep.8,2005
    idx = find(pairsCur(:,1)>pairsCur(:,2));   
    if ~isempty(idx) pairsCur(idx,:) = pairsCur(idx,[2 1]);    alignsCur(idx,:) = [single(-double(alignsCur(idx,1))) alignsCur(idx,[3 2])]; end;
    
    [common, idxOld, idxCur] = intersect(pairs,pairsCur,'rows');
    if isempty(common)
        if ~isempty(pairsCur) pairs = [pairs; pairsCur];   aligns = [aligns; alignsCur];   ratios = [ratios; load(ratioFiles{i})]; end;
    else
        ratiosCur = load(ratioFiles{i});   meansTmp{i} = load(meanFiles{i});   varsTmp{i} = load(varsFiles{i});
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

            spec1 = pairs(idxOld(j),1);   spec2 = pairs(idxOld(j),2);   updateRatio1 = meansCorrections(spec1,2)/(meansCorrections(spec1,2)+1);   updateRatio2 = meansCorrections(spec2,2)/(meansCorrections(spec2,2)+1);
            meansCorrections(spec1,:) = meansCorrections(spec1,:) + [scoresToRemove(1,1) 1];   varsCorrections(spec1) = updateRatio1*varsCorrections(spec1) + (scoresToRemove(1,1)^2)/meansCorrections(spec1,2);
            meansCorrections(spec2,:) = meansCorrections(spec2,:) + [scoresToRemove(1,2) 1];   varsCorrections(spec2) = updateRatio2*varsCorrections(spec2) + (scoresToRemove(1,2)^2)/meansCorrections(spec2,2);
        end
    end
    fprintf(1,'Finished processing %s\n',resultFiles{i});
end;

% Merge all the means/vars
proportions = zeros(size(meansTmp{1},1),numFiles);
for i=1:numFiles  proportions(:,i) = meansTmp{i}(:,2);  end;  % Get the counts
counts = sum(proportions')';   countsPos = find(counts>0);
for i=1:numFiles  proportions(countsPos,i) = proportions(countsPos,i) ./ counts(countsPos);  end;

means = zeros(size(proportions,1),1);   vars = zeros(size(proportions,1),1);
for i=1:numFiles  means = means + proportions(:,i).*meansTmp{i}(:,1);  vars = vars + proportions(:,i).*varsTmp{i}(:,1); end;

% Correct means/vars for the repeated pairs
proportions = counts(countsPos)./(counts(countsPos)-meansCorrections(countsPos,2));   
values = meansCorrections(countsPos,1)./counts(countsPos);    values_vars = varsCorrections(countsPos)./counts(countsPos);
means(countsPos) = proportions.*(means(countsPos)-values);    vars(countsPos) = proportions.*(vars(countsPos)-values_vars);

% Compute stddevs
stddevs = sqrt(vars-means.^2);
