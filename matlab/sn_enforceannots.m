function [specs, eitp, proteins] = sn_enforceannots(specs, peptides, proteins, peptidesDBS, proteinsDBS, specType, peakTol, addModMasses, minPercEI, minPercTP)
% function [specs, eitp, proteins] = sn_enforceannots(specs, peptides, proteins, peptidesDBS, proteinsDBS, specType, peakTol, addModMasses, minPercEI, minPercTP)

numSpecs = size(specs,1);   specs(:,5) = peptides;   countSame=0;
if ~isempty(peptidesDBS)
    for i=1:numSpecs 
        if ~isempty(peptidesDBS{i}) & strcmp(peptides{i},peptidesDBS{i})==1 peptidesDBS{i}=''; countSame=countSame+1; end; end;
end
fprintf(1,'  --> Got %d spectra with the same tagsearch/DBSearch peptide annotation\n',countSame);

[percs, foo, blists] = sn_getbypercs(specs, specType, peakTol, '', [], addModMasses, 0, 1, 1, 1);   eitp = percs(:,[5 1]);
[percsDBS, foo, blistsDBS] = sn_getbypercs([specs(:,1:4) peptidesDBS], specType, peakTol, '', [], addModMasses, 0, 1, 1, 1);   eitpDBS = percsDBS(:,[5 1]);
save env_annots percs blists percsDBS blistsDBS eitp eitpDBS peptides peptidesDBS;

countSwitch = 0;
for i=1:numSpecs
    if max(percsDBS(i,5),percsDBS(i,7))>max(percs(i,5),percs(i,7))  % Keep the peptide with maximum explained intensity
        specs{i,5}=peptidesDBS{i};   percs(i,:)=percsDBS(i,:);   blists(i,:)=blistsDBS(i,:);   eitp(i,:)=eitpDBS(i,:);
        peptides{i} = peptidesDBS{i};
        proteins{i} = proteinsDBS{i};
        countSwitch = countSwitch+1;
    end
    if isempty(peptides{i}) | (percs(i,5)+percs(i,7))<minPercEI | max(percs(i,1),percs(i,3))<minPercTP specs{i,5} = ''; continue; end;
    masses = sn_getmasses(peptides{i},'',[],addModMasses);   parentmass = sum(masses)+18.010564686+strcmp(specType,'msms')*1.0072763;
    masses = cumsum(masses);

    spec = specs{i,2};    specs{i,5} = peptides{i};   curBList = blists(i,:);
    if percs(i,7)>percs(i,5)  % Reverse spectrum if peaks are annotated as y
        spec(:,1) = parentmass-spec(:,1); 
        spec = sortrows(spec);
        tmpSpec = specs(i,:);   tmpSpec{1,2}=spec;
        [curPercs, foo, curBList] = sn_getbypercs(tmpSpec, specType, peakTol, '', [], addModMasses, 0, 1, 1, 1);    eitp(i,:) = curPercs([5 1]);
    end;

% fprintf(1,'%d\n',i);
% blists(i,:)
% curBList
%= blists(i,:)
    if ~isempty(masses)
%         scores = zeros(length(masses),1);   scores(curBList(:,1))=curBList(:,2)*sum(spec(:,2));
        scores = zeros(length(masses),1);   scores(curBList{1}(:,1))=spec([curBList{2}{:}],2);
        specs{i,2} = [masses' scores];
    else specs{i,2} = []; end;
end;
fprintf(1,'  --> Switched %d spectra from tagsearch to DBSearch peptide annotation (total kept: %d)\n',countSwitch,size(find(strcmp(peptides(:,1),'')==0),1));
