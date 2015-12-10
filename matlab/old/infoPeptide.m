function infoPeptide(peptide,mods,modMasses,fragsCharge)
%
%  Output an ISB like list of fragment masses for the given peptide
%

masses = getmasses3(peptide,mods,modMasses);   szMasses = size(masses,2);
bMasses = [0 cumsum(masses)+1];   yMasses = [0 cumsum(masses(szMasses:-1:1))+19];
yMasses = yMasses((szMasses+1):-1:1);
for i=1:szMasses+1
    fprintf(1,'%2d\t%6.1f\t%6.1f\t%2d\n',i-1,bMasses(i)/fragsCharge,yMasses(i)/fragsCharge,szMasses+1-i);
end
