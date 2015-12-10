function peptides = matchma_getPeptides(matchma_report)
% function peptides = matchma_getPeptides(matchma_report)
%
% Returns the best peptide match per report entry (i.e. with min # mods for maximum explained intensity)
%

numPeps = size(matchma_report,1);   peptides = cell(numPeps,1);
for i=1:numPeps
    if isempty(matchma_report{i}) continue; end;
    peptides{i} = matchma_report{i}{min(find([matchma_report{i}{:,9}]==max([matchma_report{i}{:,9}]) )),8};
end;
