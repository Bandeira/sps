function [dnstats, dninfo]= report_denovo(specs, peptides, peakTol, addModMasses)
% function dnstats = report_denovo(specs, peptides, peakTol, addModMasses)
%
%  Reports on the quality of de-novo reconstructions
%
%  dnstats(i,:) - perc. exp. int. (col.1), #FP peaks (col.2), #correct AA (col.3), #correct jumps of size 2 (col.4),
%                   #incorrect AA (col.5), #incorrect jumps of size 2 (col.6), #correct AA split by one peak (col.7)
%
%  dninfo(i,:)  - parent mass error (col.1), fragment masses offset to obtain values in dnstats(i,:) (col.2)
%
%  Problems causing unacknowledged correct sequences:
%    - Increasing peak mass error .3->.35->.4->.45->.5->.55, etc. Can be solved with shifting of mass offsets or looking at peak differences instead of masses
%
global AAmasses AAletters yProfile;
peakTolInt = round(10*peakTol);   tolRange = [-peakTolInt:peakTolInt];   szTolRange = size(tolRange,2);
jumps1 = unique(repmat(round(10*AAmasses),1,szTolRange)+repmat(tolRange,size(AAmasses,1),1));   jumps1ok = zeros(max(jumps1),1);   jumps1ok(jumps1)=1;   clear jumps1;
jumps2 = csps_getjumps(AAmasses,AAletters,2);   jumps2 = unique(repmat(round(10*jumps2),1,szTolRange)+repmat(tolRange,size(jumps2,1),1));   jumps2ok = zeros(max(jumps2),1);   jumps2ok(jumps2)=1;   clear jumps2;

if size(specs,2)==5 idxSpec=2; idxPM=3; else idxSpec=3; idxPM=5; end;
numPeps = size(peptides,1);   dnstats = zeros(numPeps,7);   dninfo = zeros(numPeps,2);
for i=1:numPeps
    if isempty(peptides{i}) | isempty(specs{i,idxSpec}) continue; end;
%     peptide = peptides{i};   spec = specs{i,idxSpec};   spec = spec(find(spec(:,1)>=57-peakTol & spec(:,1)<=specs{i,idxPM}-57+4*peakTol),:);  % Upper limit is looser to allow for parent mass errors
    peptide = peptides{i};   spec = specs{i,idxSpec};   spec = spec(find(spec(:,1)>=57-peakTol),:);

    peptideExtra='';
    if peptide(1)=='{'
        idxS=find(peptide=='{');  idxE=find(peptide=='}');   prefMass = str2num(peptide(idxS(1)+1:idxE(1)-1));   peptideExtra=peptide(idxE(1)+1:size(peptide,2));
        if prefMass>=18
            if size(idxS,2)==2 suffStr=sprintf('{%.1f}',str2num(peptide(idxS(2)+1:idxE(2)-1))+18.0); peptideExtra=peptideExtra(1:find(peptideExtra=='{')-1); else suffStr='{18.0}'; end;
            if strncmp('{18.0}',peptide,6) prefStr=''; else prefStr=sprintf('{%.1f}',prefMass-18); end;
            peptideExtra = sprintf('%s%s%s',prefStr,peptideExtra,suffStr);
        end
    end;

    masses = sn_getmasses(peptide,'',[],addModMasses);   szMasses = size(masses,2);   numPeaks = size(spec,1);   dninfo(i,1)=specs{i,idxPM}-(sum(masses)+yProfile(8,1));
    massesB = cumsum(masses(1:szMasses));   massesY = cumsum(masses(szMasses:-1:1));
    diffs = abs(repmat(massesB,numPeaks,1) - repmat(spec(:,1),1,szMasses));
%     matchesB = min(diffs)<=peakTol;    [diffs,annotsB]=min(diffs');   matchScoreB = sum(spec(find(diffs<=peakTol),2));   fpB = numPeaks - size(find(diffs<=peakTol),2);   annotsB(find(diffs>peakTol))=-1;
    if size(spec,1)>1 matchesB = min(diffs)<=peakTol; else matchesB = diffs<=peakTol; end;   [diffs,annotsB]=min(diffs');   annotsB(find(diffs>peakTol))=-1;
    [matchesB, annotsB] = extendPeakAnnots(spec, masses, peakTol, matchesB, annotsB);    % Extend matched peaks allowing for incremental mass errors
    matchScoreB = sum(spec(find(annotsB>0),2));   fpB = numPeaks - sum(matchesB);

% fprintf(1,'Offset 0: matchScoreB = %f, fpB = %d\n',matchScoreB,fpB);

    offsets = [-1 1];
    for offsetIdx=1:length(offsets)
        curMasses = cumsum(masses)+offsets(offsetIdx);
        diffs = abs(repmat(curMasses,numPeaks,1) - repmat(spec(:,1),1,szMasses));
        if size(spec,1)>1 curMatches = min(diffs)<=peakTol; else curMatches = diffs<=peakTol; end;
        [diffs,curAnnots]=min(diffs');   curAnnots(find(diffs>peakTol))=-1;
        [curMatches, curAnnots] = extendPeakAnnots(spec, masses, peakTol, curMatches, curAnnots);    % Extend matched peaks allowing for incremental mass errors

% fprintf(1,'Offset %d: matchScoreB = %f, fpB = %d\n',offsets(offsetIdx),sum(spec(find(curAnnots>0),2)),numPeaks - sum(matchesB));

        if sum(spec(find(curAnnots>0),2))>matchScoreB
            matchesB = curMatches;
            annotsB = curAnnots;
            matchScoreB = sum(spec(find(annotsB>0),2));
            fpB = numPeaks - sum(matchesB);
            dninfo(i,2) = offsets(offsetIdx);   % Fragment masses offset
        end;
    end

%     diffs = abs(repmat(massesY,numPeaks,1) - repmat(spec(:,1),1,szMasses));
%     if size(spec,1)>1 matchesY = min(diffs)<=peakTol; else matchesY = diffs<=peakTol; end;   [diffs,annotsY]=min(diffs');   matchScoreY = sum(spec(find(diffs<=peakTol),2));   fpY = numPeaks - size(find(diffs<=peakTol),2);   annotsY(find(diffs>peakTol))=-1;
%
%     if ~isempty(peptideExtra)
%         masses = sn_getmasses(peptideExtra,'',[],addModMasses);   szMasses = size(masses,2);   massesYe = cumsum(masses(szMasses:-1:1));
%         diffs = abs(repmat(massesYe,numPeaks,1) - repmat(spec(:,1),1,szMasses));
%         if size(spec,1)>1 matchesY2e = min(diffs)<=peakTol; else matchesY2e = diffs<=peakTol; end;    [diffs,annotsYe]=min(diffs');   matchScoreYe = sum(spec(find(diffs<=peakTol),2));   fpYe = numPeaks - size(find(diffs<=peakTol),2);   annotsYe(find(diffs>peakTol))=-1;
%         if matchScoreYe>matchScoreY matchesY=matchesY2e; annotsY=annotsYe; matchScoreY=matchScoreYe; fpY=fpYe; end;
%     end;
%
%     if matchScoreY>matchScoreB dnstats(i,1:2)=[matchScoreY fpY]; masses=massesY; annots=annotsY; else dnstats(i,1:2)=[matchScoreB fpB]; masses=massesB; annots=annotsB; end;
    dnstats(i,1:2)=[matchScoreB fpB]; masses=massesB; annots=annotsB;
    dnstats(i,1) = dnstats(i,1)/sum(spec(:,2));

    szAnnots = size(annots,2);
    p=1;   spec(:,1) = round(10*spec(:,1));   masses = round(10*masses);
    while p<szAnnots
        diffMass = spec(p+1,1)-spec(p,1);
        if diffMass<=size(jumps1ok,1) & jumps1ok(diffMass) massIs1aa=1; else massIs1aa=0; end;
        if diffMass<=size(jumps2ok,1) & jumps2ok(diffMass) massIs2aa=1; else massIs2aa=0; end;

        if annots(p)<0  % Check whether the current peak is incorrect
            if massIs1aa dnstats(i,5)=dnstats(i,5)+1; else if massIs2aa dnstats(i,6)=dnstats(i,6)+1; end; end;
            p=p+1; continue;
        end;

        % Look backward for the first peak
        if annots(p)==1 if p==1 dnstats(i,3)=dnstats(i,3)+1; else dnstats(i,7)=dnstats(i,7)+1; end; end
        if annots(p)==2 & p==1 dnstats(i,4)=dnstats(i,4)+1; end

        % Look to the next peak
        if annots(p+1)<0  % Check whether the next peak is incorrect
            if massIs1aa dnstats(i,5)=dnstats(i,5)+1; else if massIs2aa dnstats(i,6)=dnstats(i,6)+1; end; end;

            if (p<=szAnnots-2 & annots(p+2)>0 & annots(p+2)-annots(p)==1)
                diffMass = spec(p+2,1)-spec(p,1);   if diffMass<=size(jumps1ok,1) & jumps1ok(diffMass) dnstats(i,7)=dnstats(i,7)+1; end; % Split amino acid
            end;
            p=p+1; continue;
        end

        diff = annots(p+1)-annots(p);
        if diff==1 dnstats(i,3)=dnstats(i,3)+1; end;
        if diff==2 dnstats(i,4)=dnstats(i,4)+1; end;
        p=p+1;
    end
end

function [matches, peakAnnots] = extendPeakAnnots(spec, masses, peakTol, matches, peakAnnots);
% Extend matched peaks allowing for incremental mass errors
%
% spec    - spectrum peak masses (col.1) and scores (col.2)
% masses  - peptide masses
% matches - matched masses if no incremental mass errors are allowed
% peakAnnots - spectrum peak annotations: i-th b/y ion or -1 if not matched
%

numMasses = length(masses);   cumMasses = cumsum(masses);   numPeaks = size(spec,1);
% Extend matches from left-to-right
for i=1:numMasses-1
    if matches(i)==1 & matches(i+1)==0
        idx = find(abs(spec(:,1)-cumMasses(i))<=peakTol);
        targets = spec(idx,1)'+masses(i+1);
        if length(targets)>1 [minDiff, newIdx] = min( min( abs( repmat(spec(:,1),1,length(targets)) - repmat(targets,numPeaks,1) )' ) );
        else [minDiff, newIdx] = min( abs( repmat(spec(:,1),1,length(targets)) - repmat(targets,numPeaks,1) )' ); end;
        if minDiff<=peakTol cumMasses(i+1)=spec(newIdx,1); peakAnnots(newIdx)=i+1; matches(i+1)=1; end;
    end
end

% Extend matches from right-to-left
for i=numMasses:-1:2
    if matches(i)==1 & matches(i-1)==0
        idx = find(abs(spec(:,1)-cumMasses(i))<=peakTol);
        targets = spec(idx,1)'-masses(i);
        if length(targets)>1 [minDiff, newIdx] = min( min( abs( repmat(spec(:,1),1,length(targets)) - repmat(targets,numPeaks,1) )' ) );
        else [minDiff, newIdx] = min( abs( repmat(spec(:,1),1,length(targets)) - repmat(targets,numPeaks,1) )' ); end;
        if minDiff<=peakTol cumMasses(i-1)=spec(newIdx,1); peakAnnots(newIdx)=i-1; matches(i-1)=1; end;
    end
end

