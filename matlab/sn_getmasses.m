function masses = sn_getmasses(peptide,mods,modMasses,addModMasses)
% function masses = sn_getmasses(peptide,mods,modMasses,addModMasses)

% masses = [];
% while ~isempty(peptide)
%     if peptide(1)~='['
%         newMass = aux_getmasses(peptide(1));
%     else
%         [strMass, peptide] = strtok(peptide(2:size(peptide,2)),']');
%         newMass = str2num(strMass);  % multi-aminoacid jump
%         if isempty(newMass) 
%             newMass = aux_getmasses(strMass(1));   % equivalent amino acids [I,L] or [Q,K]
%         else
%             if abs(newMass)<=addModMasses & ~isempty(masses)
%                 masses(size(masses,2)) = masses(size(masses,2))+newMass;
%                 peptide = peptide(2:size(peptide,2));
%                 continue; 
%             end;
%         end;  
%     end;
%     peptide = peptide(2:size(peptide,2));
%     if ~isempty(peptide)
%         p=findstr(peptide(1),mods); if ~isempty(p) newMass = newMass+modMasses(p); peptide = peptide(2:size(peptide,2)); end;
%     end
%     masses = [masses newMass];
% end

masses = [];   peptide(find(peptide=='{'))='[';   peptide(find(peptide=='}'))=']';
while ~isempty(peptide)
    if peptide(1)~='['
        if peptide(1)=='(' 
            [strMass, peptide] = strtok(peptide(2:size(peptide,2)),')');
            newMass = aux_getmasses(strMass);
        else newMass = aux_getmasses(peptide(1)); end;
    else
        [strMass, peptide] = strtok(peptide(2:size(peptide,2)),']');
        idx=find(strMass==','); if ~isempty(idx) strMass=strMass(1:idx(1)-1); end;
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

function masses = aux_getmasses(peptide)
% function masses = aux_getmasses(peptide)

global AAletters AAmasses;

if isempty(peptide) masses=[]; return; end;

szPeptide = size(peptide,2);
masses = zeros(1,szPeptide);

for i=1:size(AAletters,1)
    idx = find(peptide == AAletters(i));
    masses(idx) = AAmasses(i);
end;
