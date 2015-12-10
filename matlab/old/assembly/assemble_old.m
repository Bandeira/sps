function contigs = assemble(data, aligns, mix)

numPeptides = size(data,1);

contigs = cell(numPeptides,2);
numContigs = 0;

errCount=0;

while size(aligns,1)>0
    [score, idx] = max(aligns(:,4));
    
    i = aligns(idx(1),1);
    j = aligns(idx(1),2);
    shift = aligns(idx(1),3);
    
    if mix==0 protID = 0;
    else 
        if data{i,5}==data{j,5} protID = data{i,5};
        else fprintf(1,'ERROR: Different proteins %d / %d\n',i,j);  break;
        end;
    end;
    
    posI = data{i,1};   posJ = data{j,1};   peptideI = data{i,2};   peptideJ = data{j,2};
    
    if (posI~=posJ) & (posJ-posI<0 | posJ-posI>size(peptideI,2)) 
        fprintf(1,'ERROR: (%d,%d) - posI = %d, posJ = %d, score = %.2f confidence = %2.3f, peptideI = %s, peptideJ = %s\n',i,j,posI,posJ,aligns(idx(1),4),aligns(idx(1),5),peptideI,peptideJ);
        break;
    else    
        optShift = sum(floor(getmasses(peptideI(1:(posJ-posI)))));
        if (posI~=posJ) & (shift~=optShift)
            fprintf(1,'ERROR: (%d,%d) - posI = %d, posJ = %d, score = %.2f confidence = %2.3f, peptideI = %s, peptideJ = %s, shift = %d, optShift = %d\n',i,j,posI,posJ,aligns(idx(1),4),aligns(idx(1),5),peptideI,peptideJ,shift,optShift);
%            if (posI~=posJ) errCount = errCount + 1; end;
            break;
        end;
    end;
    
    if errCount > 0 break; end;

    idxI = findContig(i,contigs);
    idxJ = findContig(j,contigs);
    
    if idxI==0 & idxJ==0 numContigs=numContigs+1; contigs{numContigs,1} = [i j];
    else if idxI==0 contigs{idxJ,1} = [contigs{idxJ,1} i]; 
        else if idxJ==0 contigs{idxI,1} = [contigs{idxI,1} j];
            else if idxI~=idxJ 
                    contigs{idxI,1} = unique([contigs{idxI,1} contigs{idxJ,1}]); 
                    contigs{idxJ,1} = contigs{numContigs,1}; 
                    numContigs=numContigs-1;
                end;
            end;
        end;
    end;
    
    
%    idx1 = find(aligns(:,1)==i);
%    idx2 = find(aligns(:,2)==j);
    idx3 = find(aligns(:,1)==j & aligns(:,2)==i);
    idxNew = setdiff([1:size(aligns,1)], idx);
%    idxNew = setdiff([1:size(aligns,1)], [idx1' idx2' idx3]);
    
%    fprintf(1,'New contig from %d to %d (%d,%d), score = %.2f confidence = %2.3f, aligns is now size %d\n',i,j,posI,posJ,aligns(idx(1),4),aligns(idx(1),5),size(aligns,1));
    
    aligns = aligns(idxNew,:);
end; % while
contigs = contigs(1:numContigs,:);

if mix==1
    for c = 1:numContigs
        pepSet = contigs{c,1};
        protSet = zeros(1,size(pepSet,2));
        for v=1:size(pepSet,2)  protSet(1,v) = data{pepSet(1,v),5};  end;
        prot = unique(protSet);
        if size(prot,2)>1 fprintf(1,'ERROR: contig set has peptides from different proteins!\n'); contigs={}; return; end;
        contigs{c,2} = prot;
    end;
end;


function idx = findContig(vertex, contigs)

idx=0;
for c=1:size(contigs,1)
    if ismember(vertex, contigs{c,1}) idx=c; return; end;
end;
