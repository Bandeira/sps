function save_bin_labels(specSet, specType, peakTol, mods, modMasses, addModMasses, filename)
% function save_bin_labels(specSet, specType, peakTol, mods, modMasses, addModMasses, filename)
%
%  Saves peak labels in binary format suitable for cpplib's LoadLabels. Each peak in specSet
%  generates one entry in filename: 0 (b-ion), 1 (y-ion), 2 (other), 3 (b _and_ y ion).
%

numSpecs=size(specSet,1);
if size(specSet,2)==5 specIdx=2; pepIdx=5; else specIdx=3; pepIdx=7; end;

numPeaks=0;   for i=1:numSpecs numPeaks = numPeaks + size(specSet{i,specIdx},1); end;
labels = int16(zeros(numPeaks,1));

baseIdx = 0;
for s=1:numSpecs
    numPeaks = size(specSet{s,specIdx},1);
    ann = annotate2(specSet{s,specIdx}, specSet{s,pepIdx}, mods, modMasses, addModMasses, strcmp(specType,'prm'), [8 11]', [6 8]', peakTol, '');
    for p=1:numPeaks
        if isempty(ann{p}) labels(baseIdx+p)=2; continue; end;
        curLabel = strcat(ann{p}{:});
        isB = ~isempty(strfind(curLabel,'b'));    isY = ~isempty(strfind(curLabel,'y'));
        if isB & isY labels(baseIdx+p)=3; continue; end;    
        if isY labels(baseIdx+p)=1; continue; end;    
        labels(baseIdx+p)=0;
    end;
    baseIdx = baseIdx + numPeaks;
end;

save_binArray(filename, labels, 'int16');
