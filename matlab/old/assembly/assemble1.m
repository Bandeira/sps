function [contigs, aligns, itersInfo] = assemble1(data, contigs, aligns, noiseProb, minConfidence, scoreFunIdx)
%
%  Assemble the contigs one at a time. At each step select the maximal scoring alignment with a 
%    confidence value greater or equal to minConfidence 

itersInfo = [];  iter=1;
while (1)
    % Filter aligns
	filtAligns = aligns(find(aligns(:,5)>=minConfidence),:);
	filtAligns = filtAligns(find(filtAligns(:,4)>100),:);
 	filtAligns = filterAligns2(contigs, filtAligns);
    filtAligns = filterCliqueSize(filtAligns, 3);
%    filtAligns = filterLimitsScore(contigs, filtAligns, 0.1);
	if (size(filtAligns,1)==0) 
        fprintf(1,'No pairwise alignments remaining after filtering\n');
        return;
	end;
	[idxOk, szWrong] = verifyAligns(contigs, data, filtAligns, 0);

    % Find highest scoring alignment
%    [maxScore, idxMax] = max(filtAligns(:,5));  % maximum confidence
    [maxScore, idxMax] = max(filtAligns(:,4));  % maximum score
    pepI = filtAligns(idxMax,1);    pepJ = filtAligns(idxMax,2);
    totalMassI = contigs{pepI,5};   totalMassJ = contigs{pepJ,5};
    shift = filtAligns(idxMax,3);
    
    fprintf(1,'Iteration %3d: Selected alignment %d (%d,%d), score = %.2f, confidence = %.2f, shift = %d\n', iter, idxMax, pepI, pepJ, filtAligns(idxMax,4), filtAligns(idxMax,5), shift);
    itersInfo = [itersInfo; size(idxOk,1) szWrong shift filtAligns(idxMax,4) filtAligns(idxMax,5)];  iter = iter+1;
    if (size(find(idxOk==idxMax),2)==0)  fprintf(1,'ORACLE: Last selected alignment is incorrect (index %d)\n', idxMax); return; end;
    
    % Merge the 2 contigs
    szContigs = size(contigs,1);
    unused = setdiff([1:szContigs], [pepI pepJ]);
    newContigs = cell(szContigs-1, size(contigs,2));
    newContigs(2:szContigs-1,:) = contigs(unused,:);

    [tmpScore, prmListOvlp] = scoreOverlap(contigs{pepI,3}(1:totalMassI,2)', contigs{pepJ,3}(1:totalMassJ,2)', shift, scoreFunIdx);
    newContigs(1,:) = mergeContigs(contigs(pepI,:), contigs(pepJ,:), shift, noiseProb, prmListOvlp, data);

    % Update aligns
    newAligns = zeros((szContigs-1)*(szContigs-2), size(aligns,2));

    newAligns = aux_selectAligns(aligns, unused', [pepI pepJ]', 1, newAligns, szContigs-1);
    newAligns = getPairAlignsC2(newContigs, [1]', [2:szContigs-1]', noiseProb, scoreFunIdx, newAligns);
    
    contigs = newContigs;
    aligns = newAligns;
end;
