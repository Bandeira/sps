function protStats = report_assembly(specs, dnSeqs, pepSeqs, components, pepMatches, pepPMatches, tags, tagPMatches, homologPMatches, topKprots, minDNTagLen, peakTol, addModMasses)
% function protStats = report_assembly(specs, dnSeqs, pepSeqs, components, pepMatches, pepPMatches, tags, tagPMatches, homologPMatches, 
%                                            topKprots, minDNTagLen, peakTol, addModMasses)
%
%  Computes statistics of assembly (coverage) and sequencing (quality)
%
%  specs - spectrum version of the reconstructed paths (as returned by masab)
%  dnSeqs - string version of the reconstructed paths (as returned by matchma)
%  pepSeqs - amino acid sequences for each entry in specs (as matched by matchma on some database)
%  components - indices of the spectra in each connected component (zero-based indices as returned by masab)
%  pepPMatches - positions where all identified spectra matched the database (as returned by findTags)
%  tagPMatches - positions where all de-novo seqs matched the database (as returned by findTags)
%  homologPMatches - positions where homologous de-novo seqs matched the database ('hand-made' in findTags output format)
%  topKprots   - reports on the top k identified proteins in pepMatches
%  minDNTagLen - only consider tags of length >= minDNTagLen
%  peakTol, addModMasses - parameters to csps_report_denovo
%
%  protStats - protein index, -1 for overall numbers (col 1)
%              coverage percentage: identified spectra, assembled spectra, matched de-novo seqs, de-novo+homologous (cols 2-5)
%              de-novo accuracy: % correct aa, % correct jumps of 2 aa, combined % (cols 6-8)
%

numSpecs = size(specs,1);
dnstats_all =csps_report_denovo(specs, pepSeqs, peakTol, addModMasses);    incrementalAApreds = zeros(1,size(dnstats_all,2));
dnlen = zeros(numSpecs,1);   cmpMatched = zeros(numSpecs,1);  
for i=1:numSpecs 
    dnlen(i) = csps_findLongestTag(dnSeqs{i},peakTol);
    cmpMatched(i)=max(pepMatches(components{i}+1))>0;    % Check whether any of the grouped spectra had a reliable annotation
end;

topKprots = min(topKprots,size(pepPMatches,1));   protsIdx = [pepPMatches{:,1}]';   protStats=zeros(topKprots+1,8);
tot_cvsums = zeros(1,5);  % Total # of aa covered by peptides(1),denovo seqs(2), homologs(3), assembled spectra(4), total proteins' length (5)
for p=1:topKprots 
    coverage = zeros(4,pepPMatches{p,3});  % coverage from single spectra (line 1), matched sequences (line 2), homologous sequences (line 3) and assembled spectra (line 4)
    coveragePeps = cell(size(coverage,2),1);  % Each position keeps a list of start/end positions for the peptides that match it. Used to count # of different overlapping peptides
    for pepIdx=1:size(pepPMatches{p,4},1) 
        idx=pepPMatches{p,4}(pepIdx,2):pepPMatches{p,4}(pepIdx,3); 
        coverage(1,idx)=coverage(1,idx)+1;
        for pos=1:size(idx,2) 
            if isempty(coveragePeps{idx(pos)}) coveragePeps{idx(pos)}=[pepPMatches{p,4}(pepIdx,2) pepPMatches{p,4}(pepIdx,3)]; continue; end; 
            if isempty(find(coveragePeps{idx(pos)}(:,1)==pepPMatches{p,4}(pepIdx,2) & coveragePeps{idx(pos)}(:,2)==pepPMatches{p,4}(pepIdx,3)))
                coveragePeps{idx(pos)}=[coveragePeps{idx(pos)}; pepPMatches{p,4}(pepIdx,2) pepPMatches{p,4}(pepIdx,3)];
            end
        end
    end;
    for pos=1:size(coverage,2) coverage(1,pos)=size(coveragePeps{pos},1); end;  % Comment this line out to get number of spectra instead of number of different overlapping peptides
    coverage(1,find(coverage(1,:)<3))=0;   
    coverage(1,:) = coverage(1,:)>0;
    protStats(p,1) = pepPMatches{p,1};

    tagPidx = find([tagPMatches{:,1}]==protsIdx(p,1));  % All tags that matched the current protein
    if ~isempty(tagPidx)
        for pepIdx=1:size(tagPMatches{tagPidx,4},1) 
            if cmpMatched(tags{tagPMatches{tagPidx,4}(pepIdx,1),1})==1 & dnlen(tags{tagPMatches{tagPidx,4}(pepIdx,1),1})>=minDNTagLen
%             if dnlen(tags{tagPMatches{tagPidx,4}(pepIdx,1),1})>=minDNTagLen
                coverage(2,tagPMatches{tagPidx,4}(pepIdx,2):tagPMatches{tagPidx,4}(pepIdx,3))=1; 
            end;
            
            % Estimate coverage of assembled identified spectra
            componentSpectra = components{tags{tagPMatches{tagPidx,4}(pepIdx,1)}}+1;
            for specIdx=1:size(componentSpectra,2) 
                if pepMatches(componentSpectra(specIdx))==protsIdx(p,1)
                    posIdx=find(pepPMatches{p,4}(:,1)==componentSpectra(specIdx)); 
                    idx=pepPMatches{p,4}(posIdx,2):pepPMatches{p,4}(posIdx,3); 
                    coverage(4,idx)=1; 
                end; 
            end;
        end;
        
        % Compute de-novo stats for all tags from components containing at least one identified spectrum and a tag of length>=minDNTagLen
        curIdx = [tags{tagPMatches{tagPidx,4}(:,1),1}]';    curIdx = curIdx(find(cmpMatched(curIdx)==1 & dnlen(curIdx)>=minDNTagLen));
%         curIdx = [tags{tagPMatches{tagPidx,4}(:,1),1}]';    curIdx = curIdx(find(dnlen(curIdx)>=minDNTagLen));
%         if size(curIdx,1)>1 curStats = sum(dnstats_all(curIdx,:)); incrementalAApreds=incrementalAApreds+curStats; end
%         if ~isempty(curIdx) protStats(p,6:8) = [curStats(1,3)/(curStats(1,3)+curStats(1,5)) curStats(1,4)/(curStats(1,4)+curStats(1,6)) (curStats(1,3)+curStats(1,4))/sum(curStats(1,3:6))]; end;
        if ~isempty(curIdx) 
            curStats = dnstats_all(curIdx,:);   if size(curIdx,1)>1 incrementalAApreds=incrementalAApreds+sum(curStats); else incrementalAApreds=incrementalAApreds+sum(curStats); end;
            protStats(p,6:8) = mean([curStats(:,3)./(curStats(:,3)+curStats(:,5)) curStats(:,4)./(curStats(:,4)+curStats(:,6)) (curStats(:,3)+curStats(:,4))./sum(curStats(:,3:6)')']); 
        end;
    end

    if ~isempty(homologPMatches)
        hmgPidx = find([homologPMatches{:,1}]==protsIdx(p,1));
        if ~isempty(hmgPidx)
            for pepIdx=1:size(homologPMatches{hmgPidx,4},1) coverage(3,homologPMatches{hmgPidx,4}(pepIdx,2):homologPMatches{hmgPidx,4}(pepIdx,3))=1; end;
        end
        coverage(3,find(coverage(2,:)==1))=0;  % Prevent double counting in homolog regions (when overlapping with well sequenced regions)
    end

    % Estimate coverage percentages
%     cvsums = sum(coverage');   protStats(p,2:5) = [cvsums(1)/pepPMatches{p,3} cvsums(4)/max(1,cvsums(1)) cvsums(2)/max(1,cvsums(1)) sum(cvsums(2:3))/max(1,cvsums(1))];    % Sequencing coverage divided by area covered by 3+ spectra
    cvsums = sum(coverage');   protStats(p,2:5) = [cvsums(1)/pepPMatches{p,3} cvsums(4)/max(1,cvsums(1)) cvsums(2)/max(1,cvsums(4)) sum(cvsums(2:3))/max(1,cvsums(4))];    % Sequencing coverage divided by area covered assembled identified spectra
    tot_cvsums = tot_cvsums + [cvsums pepPMatches{p,3}];
end
% cmpIdentified = zeros(size(cmpMatched,1),1);  cmpIdentified([tags{:,1}])=1;
% curIdx = find((cmpMatched==1 | cmpIdentified==1) & dnlen>=minDNTagLen);
curIdx = find(cmpMatched==1 & dnlen>=minDNTagLen);
if ~isempty(curIdx)
    if length(curIdx)>1 curStats = sum(dnstats_all(curIdx,:)); else curStats = dnstats_all(curIdx,:); end;

    % protStats(topKprots+1,:) = [-1 [tot_cvsums(1:2) sum(tot_cvsums(2:3))]./tot_cvsums(4) curStats(1,3)/(curStats(1,3)+curStats(1,5)) curStats(1,4)/(curStats(1,4)+curStats(1,6)) (curStats(1,3)+curStats(1,4))/sum(curStats(1,3:6))];
	% Worked when coverage had only 3 lines - protStats(topKprots+1,:) = [-1 tot_cvsums(1)/tot_cvsums(4) tot_cvsums(2)/tot_cvsums(1) sum(tot_cvsums(2:3))/tot_cvsums(1) curStats(1,3)/(curStats(1,3)+curStats(1,5)) curStats(1,4)/(curStats(1,4)+curStats(1,6)) (curStats(1,3)+curStats(1,4))/sum(curStats(1,3:6))];
	% protStats(topKprots+1,:) = [-1 tot_cvsums(1)/tot_cvsums(5) tot_cvsums(2)/tot_cvsums(1) sum(tot_cvsums(2:3))/tot_cvsums(1) curStats(1,3)/(curStats(1,3)+curStats(1,5)) curStats(1,4)/(curStats(1,4)+curStats(1,6)) (curStats(1,3)+curStats(1,4))/sum(curStats(1,3:6))];
	% protStats(topKprots+1,:) = [-1 tot_cvsums(1)/tot_cvsums(5) tot_cvsums(2)/tot_cvsums(4) sum(tot_cvsums(2:3))/tot_cvsums(4) curStats(1,3)/(curStats(1,3)+curStats(1,5)) curStats(1,4)/(curStats(1,4)+curStats(1,6)) (curStats(1,3)+curStats(1,4))/sum(curStats(1,3:6))];
	protStats(topKprots+1,:) = [-1 tot_cvsums(1)/tot_cvsums(5) tot_cvsums(4)/tot_cvsums(1) tot_cvsums(2)/tot_cvsums(4) sum(tot_cvsums(2:3))/tot_cvsums(4) curStats(1,3)/(curStats(1,3)+curStats(1,5)) curStats(1,4)/(curStats(1,4)+curStats(1,6)) (curStats(1,3)+curStats(1,4))/sum(curStats(1,3:6))];
else protStats(topKprots+1,:) = [-1 zeros(1,7)]; end;
