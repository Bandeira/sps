function [specsAnnots, resToReport, results_pathproj] = sn_load_pathproj(specsAnnots, threshExpInt, threshTP, addModMasses, ppFilename, baseFilename, components)
% function [specsAnnots, resToReport, results_pathproj] = sn_load_pathproj(specsAnnots, threshExpInt, threshTP, addModMasses, ppFilename, baseFilename, components)

if nargin<7 components=[]; specType=''; end;
if nargin<8 specType=''; end;

if size(specsAnnots,2)==5 idxPM=3; idxPep=5; else idxPM=5; idxPep=7; end;   numSpecs = size(specsAnnots,1);

htmlFilename = sprintf('report_%s_pathproj.html',baseFilename);   compsFilename = sprintf('report_%s_pathproj_components.html',baseFilename);
fid = fopen(htmlFilename,'w'); if fid<=0 fprintf(1,'ERROR opening %s!\n',htmlFilename); return; end;
fprintf(fid,'<HTML><HEAD><TITLE>pathproj search results for project %s</TITLE></HEAD><BODY>\n', baseFilename);
fprintf(fid,'<TABLE><TH>SpecIndex<TH>Annotated SpecIndex<TH>Num Mods<TH>Peptide<TH>Perc. exp. int.<TH>Perc. TP\n');

results_pathproj = sn_load_binarray(ppFilename,'int');   if isempty(results_pathproj) specsAnnots=[]; resToReport=[]; results_pathproj=[]; return; end;
results_pathproj(:,4:5) = results_pathproj(:,4:5)/10000;  % convert percentages to range [0:1]
selected = find(results_pathproj(:,4)>=threshExpInt & results_pathproj(:,5)>=threshTP);   % Selected all spectra with an eligible annotation
new = setdiff(results_pathproj(selected,2),0);  % And also all those needed to connect to a matchma/DB/original annotation
while ~isempty(setdiff(new,selected)) selected=union(selected,new);  new = setdiff(results_pathproj(selected,2),0);  end;
selected = selected(find(results_pathproj(selected,1)>0));  % Report only spectra with annotations specified by pathproj
resToReport = [selected results_pathproj(selected,:)];

[foo, idxS] = sortrows(resToReport(:,[2 1]));   resToReport = resToReport(idxS,:); % Sort by increasing level of annotation
massPos = cell(numSpecs,1);   % Indices of where masses start/end (per peptide in specsAnnots)
numResults = size(resToReport,1);
for i=1:numResults
    annotSpec = resToReport(i,3);   thisSpec = resToReport(i,1);   peptide = specsAnnots{annotSpec,idxPep};   
    
    if ~isempty(massPos{annotSpec}) curMassPos = massPos{annotSpec}; else       % Get set of mass positions for the peptide where the annotation came from
        idxStart = find((peptide>='A' & peptide<'Z') | peptide=='[')';   idxEnd = find((peptide>='A' & peptide<'Z') | peptide==']')';
        numPos = size(idxStart,1);   curMassPos = zeros(numPos,2);   curPos=1;
        for j=1:numPos   % Filter out mod masses < addModMasses from the list of start/end mass positions
            if peptide(idxStart(j))=='[' & j>1
                mass = str2num(peptide(idxStart(j)+1:idxEnd(j)-1));
                if mass<=addModMasses curMassPos(curPos-1,2) = idxEnd(j);  continue; end;  % Update previous mass entry to include this mass too
            end
            curMassPos(curPos,:) = [idxStart(j) idxEnd(j)];   curPos = curPos+1;
        end
        curMassPos = curMassPos(1:curPos-1,:);
        massPos{annotSpec} = curMassPos;
    end
    numMasses = size(curMassPos,1);
    
    if numMasses==0
        fprintf(1,'ERROR in sn_load_pathproj.m: Cannot propagate ID to %d from empty peptide %s at %d (level %d)!\n',thisSpec,peptide,annotSpec,resToReport(i,2));
        continue;
    end;
    
    % Build new peptide annotation
    modMass = specsAnnots{thisSpec,idxPM} - specsAnnots{annotSpec,idxPM};   modStr = sprintf('[%.1f]',modMass);   szModStr=size(modStr,2);
    modPos  = resToReport(i,4);   % Note that modPos is zero-based!
    if modPos==0 newPeptide = sprintf('%s%s',modStr,peptide); curMassPos = curMassPos+szModStr; curMassPos(1,1)=1; else
        if modPos==-1 newPeptide = sprintf('%s%s',peptide,modStr); curMassPos(numMasses,2) = curMassPos(numMasses,2)+szModStr; 
        else
            if modPos>=numMasses   % C-Terminal modification/extension
                newPeptide = sprintf('%s%s',peptide(1:curMassPos(numMasses,2)),modStr);
                curMassPos(numMasses,2) = curMassPos(numMasses,2)+szModStr;
            else 
                newPeptide = sprintf('%s%s%s',peptide(1:curMassPos(modPos,2)),modStr,peptide(curMassPos(modPos+1,1):curMassPos(numMasses,2))); 
                curMassPos(modPos,2) = curMassPos(modPos,2)+szModStr;   curMassPos(modPos+1:numMasses,:) = curMassPos(modPos+1:numMasses,:) + szModStr;
            end;
        end
    end
    massPos{thisSpec} = curMassPos;   specsAnnots{thisSpec,idxPep} = newPeptide;
    
    % Report results
    fprintf(fid,'<TR><TD>%.0f<TD>%.0f<TD>%.0f<TD>%s<TD>%.1f<TD>%.1f\n',thisSpec,annotSpec,resToReport(i,2),newPeptide,100*resToReport(i,5),100*resToReport(i,6));
    curEntry = annotSpec;
    while results_pathproj(curEntry,1)>=1
        fprintf(fid,'<TR><TD><TD>%.0f<TD>%.0f<TD>%s<TD>%.1f<TD>%.1f\n',curEntry,results_pathproj(curEntry,1),specsAnnots{curEntry,idxPep},100*results_pathproj(curEntry,4),100*results_pathproj(curEntry,5));
        curEntry = results_pathproj(curEntry,2);
    end;
    fprintf(fid,'<TR><TD><TD>%.0f<TD>%.0f<TD>%s<TD>%.1f<TD>%.1f\n',curEntry,results_pathproj(curEntry,1),specsAnnots{curEntry,idxPep},100*results_pathproj(curEntry,4),100*results_pathproj(curEntry,5));

end
fprintf(fid','</TABLE></BODY></HTML>');
fclose(fid);

% Output report per component
if ~isempty(components)
	fidC = fopen(compsFilename,'w'); if fidC<=0 fprintf(1,'ERROR opening %s!\n',compsFilename); return; end;
	fprintf(fidC,'<HTML><HEAD><TITLE>pathproj search results for project %s (per connected component)</TITLE></HEAD><BODY>\n', baseFilename);
    fprintf(fidC,'<p>Note: each component index is followed by number of annotated spectra / total number of spectra (per component)\n');
	fprintf(fidC,'<TABLE><TH>Index<TH>Spectrum annotation<TH>Exp Int<TH>Perc. TP<TH>Level of propagation\n');

    for cIdx=1:size(components,1)
        fprintf(fidC,'<TR><TD>Component %d<TD>%d / %d\n',cIdx,length(find(strcmp(specsAnnots(components{cIdx},idxPep),'')==0)),length(components{cIdx}));
        [foo,idxS] = sort(specsAnnots(components{cIdx},idxPep));    components{cIdx} = components{cIdx}(idxS);
        for sIdx=1:length(components{cIdx})
            curIdx = components{cIdx}(sIdx);
            if isempty(specsAnnots{curIdx,idxPep}) continue; end;
            fprintf(fidC,'<TR><TD>%d<TD>%s<TD>%.2f<TD>%.2f<TD>%d\n',curIdx,specsAnnots{curIdx,idxPep},100*results_pathproj(curIdx,4),100*results_pathproj(curIdx,5),results_pathproj(curIdx,1));
        end
    end
    
    fprintf(fidC','</TABLE></BODY></HTML>');   fclose(fidC);
end
