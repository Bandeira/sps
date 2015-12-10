function [matchStats, protStats, matchParcels] = report_csps(baseOutputFN,homglueFN, numContigs, db, refProtIdx, refProtHomologyFN, refProtSeq, peakTol, aaPerLine, contigNames)
% function [matchStats, protStats, matchParcels] = report_csps(baseOutputFN,homglueFN, numContigs, db, refProtIdx, refProtHomologyFN, refProtSeq, peakTol, aaPerLine, contigNames)
%
%  Parses homglue output and computes sequencing statistics. Mutations report is written to 'modmuts.txt' (Excel tab-separated)
%
%  matchStats   - 1 line per contig: 1-based matched prot idx (1), start/end match pos (2-3), # matched AA masses (4), # matched mass jumps (5), # mod/mut masses (6)
%  protStats    - 1 line per protein: # matched contigs (1), % coverage (2), # matched AA masses (3), # matched mass jumps (4), # mod/mut masses (5)
%  matchParcels - matched parcels per contig (csps_matchma_parcel.m output on matchma homology matches) extended to
%                  include the match-strings to the matched protein
%
%  baseOutputFN - Prefix string for output reports file names
%  homglueFN   - Name of the file containing homglue's output
%  contigs     - Contig spectra matched against the database
%  db          - Database as returned by fasta_countgis / fasta_getentries
%  refProtIdx  - If provided, indicates the index of the protein in db that is matched to the true target sequence
%  refProtHomologyFN - ClustalW homology between the true target sequence and refProtIdx (in this order)
%  refProtSeq  - True target sequence
%  peakTol     - mass tolerance used to determine mass matches/mismatches
%  aaPerLine   - Number of amino acids to report per line (default 30)
%

global AAmasses AAletters;
if nargin<5 refProtHomologyFN=''; refProtSeq=''; peakTol=0.5; aaPerLine=30; end;
if nargin<10 contigNames={}; end;
if isempty(refProtHomologyFN) refProtIdx=0; refProtCW=[]; end;

dbSize = size(db,1);
matchTxt = csps_read_csv(homglueFN,';',11);   matchTxt=matchTxt(3:size(matchTxt,1),:);
if ~isempty(refProtHomologyFN) refProtCW=csps_clustalw_load_aln(refProtHomologyFN); end;
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

    matchParcels{cIdx} = csps_matchma_parcel(matchTxt{mIdx,10});   numParcels = size(matchParcels{cIdx},1);
    matchParcels{cIdx} = [matchParcels{cIdx} cell(size(matchParcels{cIdx},1),1)];
    if protIdx~=refProtIdx
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
        matchStatsPos = [length(refProtSeq) 1];
        for pIdx=1:numParcels
            if matchParcels{cIdx}{pIdx,1}(1)=='{' continue; end;

            idx = matchcols(refProtCW(:,2),matchStats(cIdx,2)-1+[matchParcels{cIdx}{pIdx,2}:matchParcels{cIdx}{pIdx,3}]');
            idx = refProtCW(idx(find(refProtCW(idx,1)~=0)),1);
            matchStatsPos(1) = min([matchStatsPos(1) idx']);   matchStatsPos(2) = max([matchStatsPos(2) idx']);
            if isempty(idx)
                % NOTE: This code does not support whole-parcel homology-based insertions into the true sequence; partial-parcel insertions are ok (e.g. 3 AA instead of 2 will show up as a mass jump of larger mass)
                matchParcels{cIdx}{pIdx,2}=0; matchParcels{cIdx}{pIdx,3}=0;
                matchStats(cIdx,6)=matchStats(cIdx,6)+1;
                continue;
            end

            if length(idx)==1 refMass=AAmasses(find(AAletters==refProtSeq(idx))); else refMass=sum(sn_getmasses(refProtSeq(idx),'',[],0)); end;
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
                    matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} {sprintf('%c;|;',refProtSeq(idx))}];
                else
                    matchStats(cIdx,5)=matchStats(cIdx,5)+1;
                    for i=1:sz-1 matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} { sprintf('%c;.;',refProtSeq(idx(i))) }]; end;
                    matchParcels{cIdx}{pIdx,4}=[matchParcels{cIdx}{pIdx,4} {sprintf('%c;|;',refProtSeq(idx(sz)))}];
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

if refProtIdx==0
	[foo,refProtIdx] = max(protStats(:,2));   idxContigs = find(matchStats(:,1)==refProtIdx);   [foo,idxS]=sortrows(matchStats(idxContigs,2:3));   idxContigs=idxContigs(idxS);
    fprintf(1,'Selected reference protein %d: %s %s\n',refProtIdx,db{refProtIdx,1},db{refProtIdx,3});
    refProtSeq = db{refProtIdx,2};
    for i=1:length(idxContigs)
        cIdx=idxContigs(i);
        for pIdx=1:size(matchParcels{cIdx},1)
            if matchParcels{cIdx}{pIdx,1}(1)=='['
                % Output to modmuts.txt
                pos = matchStats(cIdx,2)-1+[matchParcels{cIdx}{pIdx,2} matchParcels{cIdx}{pIdx,3}];
                if pos(2)<=length(db{refProtIdx,2}) protSeq=db{refProtIdx,2}(pos(1):pos(2)); else protSeq=''; end;
                if fid>0 fprintf(fid,'%d\t%d\t%d\t%d\t%s\t%s\t%s %s\n',cIdx,refProtIdx,pos(1),pos(2),protSeq,matchParcels{cIdx}{pIdx,1},db{refProtIdx,1},db{refProtIdx,3}); end;
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

idxContigs = find(matchStats(:,1)==refProtIdx);   [foo,idxS]=sortrows(matchStats(idxContigs,2:3));   idxContigs=idxContigs(idxS);
szProtSeq = length(refProtSeq);   numBlocks = ceil(szProtSeq/aaPerLine);   blocks = cell(numBlocks,1);
blocks{1} = { ['1;;' sprintf('%c;;',refProtSeq(1:min(aaPerLine,szProtSeq)))] };
for bIdx=2:numBlocks blocks{bIdx} = { [sprintf('%d;;',(bIdx-1)*aaPerLine+1) sprintf('%c;;',refProtSeq((bIdx-1)*aaPerLine+1:min(bIdx*aaPerLine,szProtSeq)))] }; end;

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
