function masses = getmasses3(peptide,mods,modMasses)
% function masses = getmasses3(peptide,mods,modMasses)
%
%  Like getmasses but also supports [450.2] or [I,L] and modifications, e.g. mods = '*#' and modMasses = [16;32]
%

masses = [];
while ~isempty(peptide)
    if peptide(1)~='['
        newMass = getmasses(peptide(1));
    else
        [strMass, peptide] = strtok(peptide(2:size(peptide,2)),']');
        newMass = str2num(strMass);  % multi-aminoacid jump
        if isempty(newMass) newMass = getmasses(strMass(1)); end;  % equivalent amino acids
    end;
    peptide = peptide(2:size(peptide,2));
    if ~isempty(peptide)
        p=findstr(peptide(1),mods); if ~isempty(p) newMass = newMass+modMasses(p); peptide = peptide(2:size(peptide,2)); end;
    end
    masses = [masses newMass];
end