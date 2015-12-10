function masses = getMasses(peptide);
% function masses = getMasses(peptide);
%
% Get the sequence of amino acid masses for a given peptide
%
% peptide - amino acid sequence
% masses  - array of masses

global AAletters AAmasses;

if isempty(peptide) masses=[]; return; end;

szPeptide = size(peptide,2);
masses = zeros(1,szPeptide);

for i=1:size(AAletters,1)
    idx = find(peptide == AAletters(i));
    masses(idx) = AAmasses(i);
end;

