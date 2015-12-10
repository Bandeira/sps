function contigs = assemble1(data, aligns, noiseProb, bIdx, yIdx, scoreFunIdx)

contigs = [];

while (1)
    % find highest scoring alignment
    [maxScore, idx] = max(aligns(:,5));  % maximum confidence
    if maxScore<1.5 break; end;
    
    fprintf(1,'Selected alignment (%d,%d), score = %.2f, confidence = %.2f\n', aligns(idx,1), aligns(idx,2), aligns(idx,4), aligns(idx,5));
        
    % merge it and replace the 2 spectra by the new merged spectra
    ms1 = data{aligns(idx,1), 3};     totalMass1 = data{aligns(idx,1), 4};
    ms2 = data{aligns(idx,2), 3};     totalMass2 = data{aligns(idx,2), 4};
    shift = aligns(idx,3);
    
    prmIons1 = getPrmIonFrags2(totalMass1,bIdx,yIdx);
    prmIons2 = getPrmIonFrags2(totalMass2,bIdx,yIdx);

    prmScores1 = scoreprm2(ms1, [1:totalMass1]', noiseProb, 0, 1, prmIons1); % Don't care about minPeakCount for single spectra
    prmScores2 = scoreprm2(ms2, [1:totalMass2]', noiseProb, 0, 1, prmIons2);

    [tmpScore, prmList] = scoreOverlap(prmScores1, totalMass1, prmScores2, totalMass2, shift, scoreFunIdx);
    
    fprintf(1,'tmpScore = %.2f, #prms = %d\n',tmpScore, size(prmList{1,1},2));
    
    % compute the pairwise alignments between the new spectra and the existing spectra
    newPL = mergeSpec(ms1, prmIons1, ms2, prmIons2, bIdx, yIdx, shift, prmList);
    break;
end;
