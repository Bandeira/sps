function report_sps(masab_specs, matchmaData, numModStates, database, allSpecs_peptides, allSpecs_tagsMatch, allSpecs_pMatch, maxModMass, peakTol, curPrefFN, allSpecs_expint, allSpecs_bypercs, starSpecs_expint, starSpecs_bypercs, filenames, componentsFN, homologPmatches)
% function report_sps(masab_specs, matchmaData, numModStates, database, allSpecs_peptides, allSpecs_tagsMatch, allSpecs_pMatch, maxModMass, peakTol, curPrefFN,
%                              allSpecs_expint, allSpecs_bypercs, starSpecs_expint, starSpecs_bypercs, filenames, homologPmatches)
%
%
%  matchmaData - Name of the file containing the search results. However, if iscell(matchmaData) is true then
%                 this function takes the output of load_matchma from matchmaData(1,1:4).
%  numModStates - Number of modification states considered in matchmaData
%  dbSeqs - e.g. croat_db
%  allSpecs_expint,  allSpecs_bypercs  - explained intensity / by-percentages in the input spectra (to check whether input spectra were of high quality)
%  starSpecs_expint, starSpecs_bypercs - explained intensity / by-percentages in the deconvolved spectra (to check whether deconvolution was correct)
%  componentsFN    - File name of the components file output by masab (usually '../assembly/components.bla')
%  homologPmatches - Homolog regions whose resulting sequence was initially assigned to a different protein (e.g. csps_findTags.m). Same data format as other _pMatch variables.
%

if nargin<=14 componentsFN='../assembly/components.bla'; end;
if nargin<=15 homologPmatches={}; end;

components = csps_load_binListArray(componentsFN, 'int');
for cIdx = 1:size(components,1) components{cIdx}=components{cIdx}+1; end;  % Convert spec indices to one-based instead of zero-based

% Differentiate modified/unmodified peptide colors for coverage figures
save tmp_findUnmodified -v6 allSpecs_peptides;
peptides_unmod  = csps_findUnmodified(allSpecs_peptides);
idx = find([peptides_unmod{:,2}]'>0);   numPeps = size(peptides_unmod,1);
peptides_colors = [zeros(numPeps,2) ones(numPeps,1)];           % Unmodified peptides shown in blue
for i=1:length(idx) peptides_colors(idx(i),:) = [1 0 0]; end;   % Modified peptides shown in red

% Reporting results
% [report,reportScores]=csps_report_matchma(masab_specs,matchmaData, components, maxModMass, peakTol,'sps_matchma_report.html');
[report,reportScores]=csps_report_matchma(masab_specs,matchmaData, [], maxModMass, peakTol, numModStates, sprintf('%s_sps_matchma_report.html',curPrefFN));
% load croat_s2_report_sps croat_s2_report croat_s2_reportScores; report=croat_s2_report; reportScores=croat_s2_reportScores;

% tags = csps_filterLongestTag(report,3,2,maxModMass);  % Just convert almost everything to tags format
tags = csps_filterLongestTag(report,6,2,maxModMass);
if isempty(tags) fprintf(1,'(WARNING) report_sps did not complete - no contiguous matches of length 6 to the database.\n'); return; end;
[tagsMatch, protsMatch] = csps_findTags(tags(:,2),database(:,2));

% Report de-novo sequencing accuracy
pepSeqs = csps_matchma_getPeptides(report);
dnSeqs = cell(size(report,1),1); for i=1:size(report,1) dnSeqs{i}=report{i}{1,13}; end;
for cIdx = 1:size(components,1) components{cIdx}=components{cIdx}-1; end;  % Convert spec indices to zero-based instead of one-based
protStats = csps_report_assembly(masab_specs, dnSeqs, pepSeqs, components, allSpecs_tagsMatch, allSpecs_pMatch, tags, protsMatch, [], 6, 7, peakTol, maxModMass);
for cIdx = 1:size(components,1) components{cIdx}=components{cIdx}+1; end;  % Convert spec indices to one-based instead of zero-based

eval(sprintf('%s_report = report; %s_reportScores = reportScores; %s_tags = tags; %s_tagsMatch = tagsMatch; %s_protsMatch = protsMatch; %s_pepSeqs = pepSeqs; %s_dnSeqs = dnSeqs; %s_protStats = protStats;',curPrefFN,curPrefFN,curPrefFN,curPrefFN,curPrefFN,curPrefFN,curPrefFN,curPrefFN));
v=version;   v = str2num(v(1));  if v>6 v=' -V6 '; else v=''; end;
eval(sprintf('save %s_report_sps %s %s_*;',curPrefFN,v,curPrefFN));

matchedIdx = [protsMatch{:,1}]';   numMatched = size(matchedIdx,1);   reported = zeros(numMatched,1); % reported keeps track of which proteins are reported in the cycle below
for baseIdx=1:size(allSpecs_pMatch,1)
    pIdx = allSpecs_pMatch{baseIdx,1};   cIdx = find(matchedIdx==pIdx);   if isempty(cIdx) continue; end;
    reported(cIdx) = 1;

    pos = csps_getAssembledRange(allSpecs_tagsMatch, allSpecs_pMatch, [[tags{protsMatch{cIdx,4}(:,1),1}]' protsMatch{cIdx,4}(:,2:3)], pIdx, components);  % Plot assembled ranges
%     pos = protsMatch{cIdx,4}(:,2:3);  % Plot sequenced ranges

    if ~isempty(homologPmatches)
        idx = find([homologPmatches{:,1}]==pIdx);
        if ~isempty(idx) pos = [pos; homologPmatches{idx,4}(:,2:3)]; end;
    end;

    figHandle = figure('Visible','on');
    csps_plotCoverage(allSpecs_pMatch{baseIdx,3},allSpecs_pMatch{baseIdx,4}(:,2:3),1,peptides_colors(allSpecs_pMatch{baseIdx,4}(:,1),:),[],csps_mergeRanges(pos));
%     csps_plotCoverage(allSpecs_pMatch{baseIdx,3},allSpecs_pMatch{baseIdx,4}(:,2:3),1,peptides_colors(allSpecs_pMatch{baseIdx,4}(:,1),:),[],mergeContained(pos));
%     csps_plotCoverage(allSpecs_pMatch{baseIdx,3},allSpecs_pMatch{baseIdx,4}(:,2:3),1,peptides_colors(allSpecs_pMatch{baseIdx,4}(:,1),:),[],pos);
    title(sprintf('Assembled regions of %s - %s (database index %0d, ranks (%.0d,%.0d))',database{pIdx,1},database{pIdx,3},pIdx,baseIdx,cIdx),'Interpreter','none');
    xlabel('Position on the protein (in amino acids)');   ylabel('Matched MS/MS spectra');
%     saveas(figHandle,sprintf('sps_r%0d_s%0d_prot_%.0d_coverage.png',baseIdx,cIdx,pIdx),'png');
    saveas(figHandle,sprintf('%s_sps_r%0d_s%0d_prot_%.0d_coverage.fig',curPrefFN,baseIdx,cIdx,pIdx),'fig');

%     set(gcf,'Units','normalized');    set(gcf, 'Position', [0 0 1 1]);
    print('-dpng',sprintf('%s_sps_r%0d_s%0d_prot_%.0d_coverage.png',curPrefFN,baseIdx,cIdx,pIdx));
    close(figHandle);

    csps_report_coveragerange(allSpecs_pMatch{baseIdx,4},pIdx,report,components,tags,tagsMatch,curPrefFN,sprintf('%s_sps_r%0d_s%0d_prot_%.0d_report.txt',curPrefFN,baseIdx,cIdx,pIdx),allSpecs_peptides,allSpecs_expint,allSpecs_bypercs,starSpecs_expint,starSpecs_bypercs,filenames);

    if cIdx<=size(protStats,1) fprintf(1,'Protein %3d: %4.1f / %4.1f (aa seq. acc. / coverage)\n',protStats(cIdx,1),100*protStats(cIdx,6),100*protStats(cIdx,2));
    else fprintf(1,'Protein XXX: (no aa seq. acc. / coverage avalable)\n'); end;

% c=1; shifts = plotABruijnMA(gthcpt2_s1_cPAe, gthcpt2_s1_abinfo(c,:), gthcpt2_s1_report{c}{1,13}, .5); title(sprintf('gthcpt2\\_s1\\_abinfo(%0d,:), de-novo seq is %s',c,gthcpt2_s1_report{c}{1,13}));
end;

% Report proteins matched only by the contig sequences
newMatches = find(reported==0);
for newIdx=1:size(newMatches,1)
    cIdx = newMatches(newIdx);   pIdx = protsMatch{cIdx,1};   pos = protsMatch{cIdx,4}(:,2:3);  % Plot sequenced ranges

    figHandle = figure('Visible','on');
    csps_plotCoverage(protsMatch{cIdx,3},pos,1,[],[],[]);
    title(sprintf('Ranks (none,%.0d), database index %0d : %s - %s',cIdx,pIdx,database{pIdx,1},database{pIdx,3}),'Interpreter','none');
    saveas(figHandle,sprintf('%s_sps_rXX_s%0d_prot_%.0d_coverage.fig',curPrefFN,cIdx,pIdx),'fig');
%     print('-djpeg',sprintf('sps_rXX_s%0d_prot_%.0d_coverage.jpg',cIdx,pIdx));
    print('-dpng',sprintf('%s_sps_rXX_s%0d_prot_%.0d_coverage.png',curPrefFN,cIdx,pIdx));
    close(figHandle);

    specIdx = [components{[tags{protsMatch{cIdx,4}(:,1),1}]'}]';
    specIdx = specIdx(find(allSpecs_tagsMatch(specIdx)==0));  % Fake matches to this protein for all spectra that did not match anything. Helps distiguish from contigs assembling spectra that matched other proteins.
    csps_report_coveragerange(specIdx,pIdx,report,components,tags,tagsMatch,curPrefFN,sprintf('%s_sps_rXX_s%0d_prot_%.0d_report.txt',curPrefFN,cIdx,pIdx),allSpecs_peptides,allSpecs_expint,allSpecs_bypercs,starSpecs_expint,starSpecs_bypercs, filenames);
    if cIdx<=size(protStats,1) fprintf(1,'Protein %3d: %4.1f / %4.1f (aa seq. acc. / coverage)\n',protStats(cIdx,1),100*protStats(cIdx,6),100*protStats(cIdx,2));
    else fprintf(1,'Protein XXX: (no aa seq. acc. / coverage avalable)\n'); end;
end;

% Report contig sequences that did not match any protein
orphan = find(tagsMatch==0);   szOrphan = size(orphan,1);   fid = fopen(sprintf('%s_sps_orphan_report.txt',curPrefFN),'w');
if fid<=0 fprintf(1,'Error opening sps_orphan_report.txt!\n'); else
    fprintf(fid,'%0d contig sequences did not match any protein.\n\n',szOrphan);
    for oIdx=1:szOrphan
        compIdx = orphan(oIdx);   numSpecs = size(components{compIdx},2);
        fprintf(fid,'Component %.0f\n',compIdx);
        fprintf(fid,'   De-novo sequence    : %s \n',report{compIdx}{1,13});
        fprintf(fid,'   # assembled spectra : %.0f \n',numSpecs);

        [foo, idxS] = sort(allSpecs_peptides(components{compIdx}));
        for specIdx=1:numSpecs fprintf(fid,'   (protein %4.0d) %s\n',allSpecs_tagsMatch(components{compIdx}(idxS(specIdx))),allSpecs_peptides{components{compIdx}(idxS(specIdx))}); end;
        fprintf(fid,'\n');
    end;
    fclose(fid);
end;

