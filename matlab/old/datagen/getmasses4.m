function masses = getmasses4(peptide,mods,modMasses,addModMasses)
% function masses = getmasses4(peptide,mods,modMasses,addModMasses)
%
%  Like getmasses3 (supports [450.2] or [I,L] and modifications, e.g. mods = '*#' and modMasses = [16;32])
%    but when addModMasses>0 then T[mod] generates only one mass T+mod whenever abs(mod)<=addModMasses
%

masses = [];
while ~isempty(peptide)
    if peptide(1)~='['
        newMass = getmasses(peptide(1));
    else
        [strMass, peptide] = strtok(peptide(2:size(peptide,2)),']');
        newMass = str2num(strMass);  % multi-aminoacid jump
        if isempty(newMass) 
            newMass = getmasses(strMass(1));   % equivalent amino acids [I,L] or [Q,K]
        else
            if abs(newMass)<=addModMasses & ~isempty(masses)
                masses(size(masses,2)) = masses(size(masses,2))+newMass;
                peptide = peptide(2:size(peptide,2));
                continue; 
            end;
        end;  
    end;
    peptide = peptide(2:size(peptide,2));
    if ~isempty(peptide)
        p=findstr(peptide(1),mods); if ~isempty(p) newMass = newMass+modMasses(p); peptide = peptide(2:size(peptide,2)); end;
    end
    masses = [masses newMass];
end
