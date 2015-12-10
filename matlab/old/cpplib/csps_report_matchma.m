function [report, matchScores, results] = report_matchma(specs, matchmaData, vSets, addModMasses, peakTol, numModStates, htmlFN)
% function [report, matchScores, results] = report_matchma(specs, matchmaData, vSets, addModMasses, peakTol, numModStates, htmlFN)
%
%  Creates an HTML report showing the results of running matchma on specs/vSets
%
%  matchmaData - Name of the file containing the search results. However, if iscell(matchmaData) is true then
%                 this function takes the output of csps_load_matchma from matchmaData(1,1:4).
%  vSets - sets of indices of assembled spectra per component (only used for reporting component size)
%
%  report - each entry contains the same informatation that is output into htmlFN
%  matchScores - each entry contains [% exp. score, %TP on DB peptide]
%  results = as returned by csps_load_matchma
%

global AAmasses AAletters;
fid=fopen(htmlFN,'w');  if(fid<=0) fprintf(1,'Error opening %s!\n',htmlFN); return; end;
if size(specs,2)==5 idxSpec=2; idxPM=3; idxPep=5; else idxSpec=3; idxPM=5; idxPep=7; end;
if isempty(vSets) vSets=cell(size(specs,1),1); for i=1:size(specs,1) vSets{i}=i; end; end;

if iscell(matchmaData)  % Matchma results already loaded
	peptidesDB    = matchmaData{1,1};
    peptidesMatch = matchmaData{1,2};
    peptidesPos   = matchmaData{1,3};
    results       = matchmaData{1,4};
	numEntries = size(peptidesDB,1);
else
	numEntries = size(vSets,1);
	[peptidesDB, peptidesMatch, peptidesPos, results] = csps_load_matchma(matchmaData, ';', numEntries, numModStates);
%fprintf(1,'105: (%d,%d), (%d,%d), numEntries = %d\n',size(results{105,1}),size(results{105,2}),numEntries);
	numEntries = size(peptidesDB,1); if numEntries~=size(vSets,1) fprintf(1,'Error: incorrect number of entries in %s!\n',matchmaData); fclose(fid); return; end;
	fprintf(1,'Finished reading %s (%.0f entries)\n',matchmaData,numEntries);
end;

fprintf(fid,'<HTML><HEAD><TITLE>Matchma search results</TITLE></HEAD><BODY>\n');
fprintf(fid,'<TABLE><TH>Count<TH>SpecIndex<TH>Num Specs<TH>Num Mods<TH>b/y match<TH>Len<TH>Peptide<TH>Peptide match<TH>Perc. exp. int.<TH>Perc. TP<TH>TP<TH>FP<TH>Denovo peptide\n');
report=cell(numEntries,1);   reportStrings = cell(numEntries,1);   entryScore = zeros(numEntries,1);   matchScores = zeros(numEntries,2);
for i=1:numEntries
    prefix = sprintf('<TD>%d<TD>%d',i,size(vSets{i},2));
    prefixRep = {i,size(vSets{i},2)};
    denovoPep = csps_denovoPeptide(specs{i,idxSpec},1,AAmasses,AAletters,min(peakTol,0.47));  %  0.47 max peakTol avoids Q/K/E confusions
    for j=1:size(peptidesDB,2)
        curSpec = specs(i,:);
        if isempty(results{i,j})
            percScore=0; numTP=0; percTP=0; numFP=size(curSpec{1,idxSpec},1); dbLen=0; matchType='';
        else
            curSpec{1,idxPep} = match2peptide(peptidesMatch{i,j});   matchLen = size(sn_getmasses(curSpec{1,idxPep},'',[],addModMasses),2);  % matchLen = # of matched prefix masses
            matchType=results{i,j}{1,1};
            if matchType=='y'
                curSpec{1,idxSpec}(:,1) = curSpec{1,idxSpec}(:,1)+18;
                if strncmp(curSpec{1,idxPep},'[18.0]',6)==1 curSpec{1,idxPep}=curSpec{1,idxPep}(7:size(curSpec{1,idxPep},2)); end;
                curSpec = repmat(curSpec,2,1);   curSpec{1,idxSpec} = specs{i,idxSpec};
            end;
            idx = find(curSpec{1,idxSpec}(:,1)>50 & curSpec{1,idxSpec}(:,1)<curSpec{1,idxPM}-50);  % remove trailing endpoints
            curSpec{1,idxSpec} = curSpec{1,idxSpec}(idx,:);   if size(curSpec,1)>1 curSpec{2,idxSpec} = curSpec{2,idxSpec}(idx,:); end;
            numPeaks = size(curSpec{1,idxSpec},1);

            [expInt, expPeaks] = csps_getsnr2(curSpec,'prm','',[],addModMasses,peakTol,[]);
            bypercs = sn_getbypercs(curSpec,'prm',peakTol,'',[],addModMasses,0,0);
            if size(curSpec,1)>1 bestYmatch = find(bypercs(:,7)==max(bypercs(:,7))); expPeaks=expPeaks(bestYmatch,:);   bypercs=bypercs(bestYmatch,:);  end;
            pPos=find(peptidesDB{i,j}=='.'); dbLen = max(pPos)-min(pPos)-1;  % dbLen = # amino acids in database sequence
            if results{i,j}{1,1}=='b' percScore=bypercs(1,5); percTP=bypercs(1,1); else percScore=bypercs(1,7); percTP=bypercs(1,3); end;
            numTP  = percTP*bypercs(1,9);
            numFP  = expPeaks(1,7)*expPeaks(1,8);
        end;
        reportStrings{i} = [reportStrings{i};{sprintf('%s<TD>%d<TD>%s<TD>%d<TD align=left nowrap>%s<TD align=left nowrap>%s<TD>%4.1f<TD>%4.1f<TD>%.0f<TD>%.0f<TD align=left nowrap>%s\n',prefix,j-1,matchType,dbLen,peptidesDB{i,j},peptidesMatch{i,j},100*percScore,100*percTP,numTP,numFP,denovoPep)}];
        report{i} = [report{i}; [prefixRep {prefix,j-1,matchType,dbLen,peptidesDB{i,j},peptidesMatch{i,j},100*percScore,100*percTP,numTP,numFP,denovoPep}]];
        oldValue = entryScore(i);
        entryScore(i) = max(entryScore(i),percScore+percTP);
        if entryScore(i)>sum(matchScores(i,:)) matchScores(i,:) = [percScore percTP]; end;
        if entryScore(i)>=2-0.0001 break; end;
        prefix = '<TD><TD>';   denovoPep = '';
    end
end

[foo, idxS] = sort(entryScore);   idxS = idxS(numEntries:-1:1);  reportStrings=reportStrings(idxS,:);  % sort by decreasing entryScore
for i=1:numEntries
    fprintf(fid,'<TR align=CENTER><TD>%d',i); fprintf(fid,reportStrings{i}{1});
    for j=2:size(reportStrings{i},1) fprintf(fid,'<TR align=CENTER><TD>'); fprintf(fid,reportStrings{i}{j}); end;
    fprintf(fid,'<TR><TD height=5>\n');
end;

fprintf(fid','</TABLE></BODY></HTML>');
fclose(fid);

function strPep = match2peptide(strMatch)
%
%  Converts a matchma peptide match string to a peptide annotation string (i.e. suitable for csps_getsnr2)
%  {} -> [], (ABC) -> ABC, [A,B,C] -> [A]

if isempty(strMatch) strPep=strMatch; return; end;
% remove ()
strPep = strMatch(find(strMatch~='(' & strMatch~=')'));

% [A,B,C] -> [A]
posS = find(strPep=='[')';   posE = find(strPep==']')';
for i=1:size(posS,1)
%    idx=min(find(strPep(posS(i)+1:posE(i)-1)==','));
%    strPep(posS(i)+idx:posE(i)-1) = ' ';

    % [A,B,C] -> B[C]
    idx=find(strPep(posS(i)+1:posE(i)-1)==',');
    strPep(posS(i):posS(i)+min(idx))=' ';   strPep(posS(i)+max(idx))='[';
end
strPep = strPep(find(strPep~=' '));

% {} -> []
strPep(find(strPep=='{'))='[';   strPep(find(strPep=='}'))=']';

