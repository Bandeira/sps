function [matchStats, protStats, matchParcels] = csps_report(baseOutputFN,homglueFN, numContigs, db, refProtIdx, trueProtIdx, trueProtHomologyFN, trueProtSeq, peakTol, aaPerLine, contigNames)
% function [matchStats, protStats, matchParcels] = csps_report(baseOutputFN,homglueFN, numContigs, db, refProtIdx, trueProtIdx, trueProtHomologyFN, trueProtSeq, peakTol, aaPerLine, contigNames)
%
%  Parses homglue output and computes sequencing statistics. Mutations report is written to 'modmuts.txt' (Excel tab-separated)
%
%  matchStats   - 1 line per contig: 1-based matched prot idx (1), start/end match pos (2-3), # matched AA masses (4), # matched mass jumps (5), # mod/mut masses (6)
%  protStats    - 1 line per protein: # matched contigs (1), % coverage (2), # matched AA masses (3), # matched mass jumps (4), # mod/mut masses (5)
%  matchParcels - matched parcels per contig (matchma_parcel.m output on matchma homology matches) extended to
%                  include the match-strings to the matched protein
%
%  baseOutputFN - Prefix string for output reports file names
%  homglueFN   - Name of the file containing homglue's output
%  contigs     - Contig spectra matched against the database
%  db          - Database as returned by fasta_countgis / fasta_getentries
%  trueProtIdx  - If provided, indicates the index of the protein in db that is matched to the true target sequence
%  trueProtHomologyFN - ClustalW homology between the true target sequence and trueProtIdx (in this order)
%  trueProtSeq  - True target sequence
%  peakTol     - mass tolerance used to determine mass matches/mismatches
%  aaPerLine   - Number of amino acids to report per line (default 30)
%

global AAmasses AAletters;
if nargin<6 trueProtHomologyFN=''; trueProtSeq=''; peakTol=0.5; aaPerLine=30; end;
if nargin<11 contigNames={}; end;
if isempty(trueProtHomologyFN) trueProtIdx=-trueProtIdx; trueProtCW=[]; end;

dbSize = size(db,1);
matchTxt = csps_load_csv(homglueFN,';',11);   matchTxt=matchTxt(2:size(matchTxt,1),:);
if ~isempty(trueProtHomologyFN) trueProtCW=csps_load_clustalw_aln(trueProtHomologyFN); end;
fid = fopen(strcat(baseOutputFN,'_modmuts.txt'),'w');   if fid>0 fprintf(fid,'ContigIdx\tProteinIdx\tStart aa\tEnd aa\tMatched\tMassDiff\tProteinID\n'); end;

matchStats = zeros(numContigs,6);   
matchParcels = cell(numContigs,1);
protStats = zeros(dbSize,5);        
protsCvg = cell(dbSize,1);   for pIdx=1:dbSize protsCvg{pIdx}=zeros(1,length(db{pIdx,2})); end;
for mIdx=1:size(matchTxt,1)
    cIdx = str2num(matchTxt{mIdx,2})+1;
    matchStats(cIdx,1) = str2num(matchTxt{mIdx,6})+1;   protIdx = matchStats(cIdx,1);
    matchStats(cIdx,2) = str2num(matchTxt{mIdx,7})+1;
    matchStats(cIdx,3) = str2num(matchTxt{mIdx,8});
    protStats(protIdx,1) = protStats(protIdx,1)+1;
    
    matchParcels{cIdx} = aux_matchma_parcel(matchTxt{mIdx,11});   numParcels = size(matchParcels{cIdx},1);
    matchParcels{cIdx} = [matchParcels{cIdx} cell(size(matchParcels{cIdx},1),1)];
    if protIdx~=trueProtIdx
        for pIdx=1:numParcels
            if matchParcels{cIdx}{pIdx,1}(1)=='{' continue; end;
            
            if matchParcels{cIdx}{pIdx,1}(1)=='(' 
                matchStats(cIdx,5)=matchStats(cIdx,5)+1;   sz = length(matchParcels{cIdx}{pIdx,1});
                for i=2:sz-2 matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} { sprintf('%c;.;',matchParcels{cIdx}{pIdx}(i)) }]; end;
                matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} {sprintf('%c;|;',matchParcels{cIdx}{pIdx}(sz-1))}];
                continue; 
            end;
            
            if matchParcels{cIdx}{pIdx,1}(1)=='[' 
                matchStats(cIdx,6)=matchStats(cIdx,6)+1;   idxC = find(matchParcels{cIdx}{pIdx,1}==',');   sz = idxC(2)-idxC(1)-1;

                if sz==1 matchParcels{cIdx}{pIdx,4}= { sprintf('%s;|;',matchParcels{cIdx}{pIdx}(2:idxC(1)-1)) } ;
                else 
                    matchParcels{cIdx}{pIdx,4}= { sprintf('%s;;',matchParcels{cIdx}{pIdx}(2:idxC(1)-1)) } ;
                    matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} repmat({';;'},1,sz-2)];
                    matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} { ';|;' }];
                end;
                continue; 
            end;
            
            matchStats(cIdx,4)=matchStats(cIdx,4)+1;
            matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} {sprintf('%s;|;',matchParcels{cIdx}{pIdx})}];
        end
    else
        parcelPosOffset=0;   % offset between parcel indices and protein-match indices; non-zero when protein homology contains aa insertions/deletions
        matchStatsPos = [length(trueProtSeq) 1];
        for pIdx=1:numParcels
            if matchParcels{cIdx}{pIdx,1}(1)=='{' continue; end;
            
            idx = aux_matchCols(trueProtCW(:,2),matchStats(cIdx,2)-1+[matchParcels{cIdx}{pIdx,2}:matchParcels{cIdx}{pIdx,3}]');
            idx = trueProtCW(idx(find(trueProtCW(idx,1)~=0)),1);
            matchStatsPos(1) = min([matchStatsPos(1) idx']);   matchStatsPos(2) = max([matchStatsPos(2) idx']);
            if isempty(idx)
                % NOTE: This code does not support whole-parcel homology-based insertions into the true sequence; partial-parcel insertions are ok (e.g. 3 AA instead of 2 will show up as a mass jump of larger mass)
                matchParcels{cIdx}{pIdx,2}=0; matchParcels{cIdx}{pIdx,3}=0; 
                matchStats(cIdx,6)=matchStats(cIdx,6)+1; 
                continue; 
            end

            if length(idx)==1 refMass=AAmasses(find(AAletters==trueProtSeq(idx))); else refMass=sum(sn_getmasses(trueProtSeq(idx),'',[],0)); end;
            parcelPosOffset = parcelPosOffset+length(idx)-(matchParcels{cIdx}{pIdx,3}-matchParcels{cIdx}{pIdx,2}+1);
            sameLen = (length(idx)==(matchParcels{cIdx}{pIdx,3}-matchParcels{cIdx}{pIdx,2}+1));
            matchParcels{cIdx}{pIdx,2}=matchParcels{cIdx}{pIdx,2}+parcelPosOffset;
            matchParcels{cIdx}{pIdx,3}=matchParcels{cIdx}{pIdx,2}+length(idx)-1;
            
            if sameLen==1 & length(matchParcels{cIdx}{pIdx,1})==1
                curMass=AAmasses(find(AAletters==matchParcels{cIdx}{pIdx,1}));
                if abs(curMass-refMass)<=peakTol matchStats(cIdx,4)=matchStats(cIdx,4)+1; else matchStats(cIdx,6)=matchStats(cIdx,6)+1; end;
                matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} {sprintf('%s;|;',matchParcels{cIdx}{pIdx,1})}];
                continue;
            end
            
            % Process ()/[] parcels
            curMass = sum(sn_getmasses(matchParcels{cIdx}{pIdx,1},'',[],0));   sz = length(idx);
            if abs(curMass-refMass)<=peakTol 
                if matchParcels{cIdx}{pIdx,1}(1)=='[' & matchParcels{cIdx}{pIdx,2}==matchParcels{cIdx}{pIdx,3}
                    matchStats(cIdx,4)=matchStats(cIdx,4)+1;
                    matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} {sprintf('%c;|;',trueProtSeq(idx))}];
                else 
                    matchStats(cIdx,5)=matchStats(cIdx,5)+1;
                    for i=1:sz-1 matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} { sprintf('%c;.;',trueProtSeq(idx(i))) }]; end;
                    matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} {sprintf('%c;|;',trueProtSeq(idx(sz)))}];
                end;
            else 
                matchStats(cIdx,6)=matchStats(cIdx,6)+1;   idxC = find(matchParcels{cIdx}{pIdx,1}==',');
                if sz==1 matchParcels{cIdx}{pIdx,4}= { sprintf('%s;|;',matchParcels{cIdx}{pIdx}(2:idxC(1)-1)) } ;
                else
                    if matchParcels{cIdx}{pIdx,1}(1)=='(' 
                        if length(matchParcels{cIdx}{pIdx,1})==sz+2 txt=matchParcels{cIdx}{pIdx,1}(2:sz+1); else txt=repmat('X',1,sz); end;
                        for i=1:sz-1 matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} {sprintf('%c;.;',txt(i))}]; end;
                        matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} {sprintf('%c;|;',txt(sz))}];
                    else
                        matchParcels{cIdx}{pIdx,4}= { sprintf('%s;;',matchParcels{cIdx}{pIdx}(2:idxC(1)-1)) } ;
                        matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} repmat({';;'},1,sz-2)];
                        matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} { ';|;' }];
                    end;
                end;
            end;
        end
        matchStats(cIdx,2:3) = matchStatsPos;
    end
    protsCvg{protIdx}(matchStats(cIdx,2):matchStats(cIdx,3)) = 1;  % Set of covered amino acids
    protStats(protIdx,3:5) = protStats(protIdx,3:5) + matchStats(cIdx,4:6);
end

for protIdx=1:dbSize
    protStats(protIdx,2) = length(find(protsCvg{protIdx}==1)) / length(protsCvg{protIdx});
end

if isempty(trueProtHomologyFN) trueProtIdx=-trueProtIdx; end;
if isempty(trueProtHomologyFN) | trueProtIdx==0
%     if trueProtIdx==0
% %       [foo,trueProtIdx] = max(protStats(:,2));
% %       [foo,trueProtIdx] = max(protStats(:,1));
%         [foo,trueProtIdx] = max(sum(protStats(:,3:4)')');
%     end;
    trueProtIdx = refProtIdx;
    idxContigs = find(matchStats(:,1)==trueProtIdx);   [foo,idxS]=sortrows(matchStats(idxContigs,2:3));   idxContigs=idxContigs(idxS);
    fprintf(1,'Selected reference protein %d: %s %s\n',trueProtIdx,db{trueProtIdx,1},db{trueProtIdx,3});
    trueProtSeq = db{trueProtIdx,2};
    for i=1:length(idxContigs)
        cIdx=idxContigs(i);
        for pIdx=1:size(matchParcels{cIdx},1)
            if matchParcels{cIdx}{pIdx,1}(1)=='['
                % Output to modmuts.txt
                pos = matchStats(cIdx,2)-1+[matchParcels{cIdx}{pIdx,2} matchParcels{cIdx}{pIdx,3}];
                if pos(2)<=length(db{trueProtIdx,2}) protSeq=db{trueProtIdx,2}(pos(1):pos(2)); else protSeq=''; end;
                if fid>0 fprintf(fid,'%d\t%d\t%d\t%d\t%s\t%s\t%s %s\n',cIdx,trueProtIdx,pos(1),pos(2),protSeq,matchParcels{cIdx}{pIdx,1},db{trueProtIdx,1},db{trueProtIdx,3}); end;
            end;
        end
    end
end
if fid>0 fclose(fid); end;

%
%  Output matches to semicolon-separated text file
%

fidOvlp = fopen(strcat(baseOutputFN,'_homglue_overlaps.txt'),'w');    if fidOvlp<=0 fprintf(1,'ERROR opening %s_homglue_overlaps.txt!\n',baseOutputFN); return; end;
% fidOvlp = 1;

idxContigs = find(matchStats(:,1)==trueProtIdx);   [foo,idxS]=sortrows(matchStats(idxContigs,2:3));   idxContigs=idxContigs(idxS);
szProtSeq = length(trueProtSeq);   numBlocks = ceil(szProtSeq/aaPerLine);   blocks = cell(numBlocks,1);
blocks{1} = { ['1;;' sprintf('%c;;',trueProtSeq(1:min(aaPerLine,szProtSeq)))] };
for bIdx=2:numBlocks blocks{bIdx} = { [sprintf('%d;;',(bIdx-1)*aaPerLine+1) sprintf('%c;;',trueProtSeq((bIdx-1)*aaPerLine+1:min(bIdx*aaPerLine,szProtSeq)))] }; end;

for i=1:length(idxContigs)
    cIdx=idxContigs(i);    bIdx = floor((matchStats(cIdx,2)-1)/aaPerLine)+1;   outPos = rem(matchStats(cIdx,2)-1,aaPerLine)+1;
    if isempty(contigNames) cName=sprintf('%d',cIdx); else cName=sprintf('%s',contigNames{cIdx}); end;
    curLine = strcat(cName,sprintf('%s;|;',repmat(';;',1,rem(matchStats(cIdx,2)-1,aaPerLine))));
%     if matchStats(cIdx,2)>1 & outPos==1 blocks{bIdx} = [blocks{bIdx}; {curLine}]; bIdx=bIdx+1; curLine = sprintf('%d;;',cIdx); end;
    
    lastBreak='';
    for pIdx=1:size(matchParcels{cIdx},1)
        for i=1:size(matchParcels{cIdx}{pIdx,4},2)
            curLine=strcat(curLine,matchParcels{cIdx}{pIdx,4}{i});
            outPos=outPos+1;
            lastBreak = matchParcels{cIdx}{pIdx,4}{i};   lastBreak = lastBreak(length(lastBreak)-1);   if lastBreak==';' lastBreak=''; end;
            if outPos>aaPerLine outPos=1; blocks{bIdx} = [blocks{bIdx}; {curLine}]; bIdx=bIdx+1; curLine = sprintf('%s;%s;',cName,lastBreak); end;
        end
    end
    if ~strcmp(curLine,sprintf('%s;%c;',cName,lastBreak)) blocks{bIdx}=[blocks{bIdx}; {curLine}]; end;
end

% Output blocks
for bIdx=1:numBlocks
    for lineIdx=1:size(blocks{bIdx},1)
        fprintf(fidOvlp,'%s\n',blocks{bIdx}{lineIdx});
    end
    fprintf(fidOvlp,'\n');
end
fclose(fidOvlp); 

function parcels = aux_matchma_parcel(str)
% function parcels = aux_matchma_parcel(str)
%
%  Divides a string into {}, (), [] blocks and individual matched amino acids.
%
%  parcels - as many lines as there are parcels, each line is a parcel string (1), AA match start/end pos (2-3)
%

idxOpen = find(str=='{' | str=='(' | str=='[');
idxClose = find(str=='}' | str==')' | str==']');
numBlocks=length(idxOpen);   blocksLen = sum(idxClose-idxOpen+1);   numParcels = numBlocks+length(str)-blocksLen;

parcels = cell(numParcels,3);   curPos=1;   pIdx=1;   curAApos=1;
for bIdx=1:numBlocks
    if idxOpen(bIdx)>curPos
        numAA = idxOpen(bIdx)-curPos;
        for i=1:numAA
            parcels{pIdx,1}=str(curPos+i-1);   
            parcels{pIdx,2}=curAApos;   parcels{pIdx,3}=curAApos;   curAApos=curAApos+1;
            pIdx = pIdx+1;
        end
    end
    parcels{pIdx,1}=str(idxOpen(bIdx):idxClose(bIdx));
    if parcels{pIdx,1}(1)~='{' 
        parcels{pIdx,2}=curAApos;
        if parcels{pIdx,1}(1)=='(' parcels{pIdx,3}=curAApos+length(parcels{pIdx,1})-3; else
            idx = find(parcels{pIdx,1}==',');   txt = parcels{pIdx,1}(idx(1)+1:idx(2)-1);
            parcels{pIdx,3}=curAApos+length(txt)-1;
        end
        curAApos=parcels{pIdx,3}+1;
    end;
    pIdx = pIdx+1;   curPos = idxClose(bIdx)+1;
end

if isempty(idxOpen) | idxClose(numBlocks)<length(str)
    numAA = length(str)-curPos+1;
    for i=1:numAA
        parcels{pIdx}=str(curPos+i-1);
        parcels{pIdx,2}=curAApos;   parcels{pIdx,3}=curAApos;   curAApos=curAApos+1;
        pIdx = pIdx+1;
    end
end

function idx=aux_matchCols(c1,c2)
% function idx=matchCols(c1,c2)
%
%  idx is the index in c1 of elements in c2
%

matches = zeros(size(c1,1),1);
for i=1:size(c2,1)
    matches(find(c1==c2(i)))=1;
end
idx=find(matches==1);
