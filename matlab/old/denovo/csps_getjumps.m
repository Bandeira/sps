function [masses, letters] = getjumps(inMasses, inLetters, k)
% function [masses, letters] = getjumps(inMasses, inLetters, k)
%
%  Like get jumps (builds a structure with all the jumps of k masses) but letters is
%   now a cell variable with a list of amino acid combinations that yield the considered mass.
%

if k==0 masses=[]; letters=''; return; end;

baseMasses = inMasses;   resolution = .1;
useMasses = round(inMasses/resolution);    masses = unique(useMasses);    letters = cell(size(masses));   
for iter=1:k
	for i=1:size(useMasses,1)
        idx = find(masses==useMasses(i));
        if iter<2
            if isempty(letters{idx,1})  letters{idx,1} = inLetters(i);
            else
                letters{idx,1} = sprintf('[%s,%s]', letters{idx,1}, inLetters(i));
            end
        else
            letters{idx,1} = sprintf('[%.1f]',masses(idx)/10);
        end;
    end;
    if iter<k
        inMasses = unique(reshape( repmat(inMasses, 1, size(baseMasses,1)) + repmat(baseMasses', size(inMasses,1), 1) , size(inMasses,1)*size(baseMasses,1), 1));
        useMasses = setdiff(unique(round(inMasses/resolution)), masses);
        masses = [masses; useMasses];
    end;
end;
masses = masses * resolution;
