function compsReported = report_coveragerange(specIndices, curProtIdx, matchmaReport, components, tags, tagsMatch, curFNpref, reportFN, allPeptides, ctgs_expint, ctgs_bypercs, cPAe_expint, cPAe_bypercs, filenames)
% function compsReported = report_coveragerange(specIndices, curProtIdx, matchmaReport, components, tags, tagsMatch, curFNpref, reportFN, allPeptides, 
%                                        ctgs_expint, ctgs_bypercs, cPAe_expint, cPAe_bypercs, filenames)
%
%  Reports alignment/grouping/assembly results for a set of spectra.
%
%  specIndices - Indices of spectra to report about
%  curProtIdx  - protein where these spectra are supposed to come from (protein index in the database used for tags/tagsMatch).
%  matchmaReport  - Used to obtain identified sequences for components with spectra in specIndices
%  tags/tagsMatch - Location of all individual peptides on the target proteins (as returned by findTags)
%  allPeptides    - Peptide annotations for every spectrum
%  ctgs_expint, ctgs_bypercs - explained intensity / by-percentages in the input spectra (to check whether input spectra were of high quality)
%  cPAe_expint, cPAe_bypercs - explained intensity / by-percentages in the deconvolved spectra (to check whether deconvolution was correct)
%
%  compsReported - indices of all components reported

% components = load_binListArray('../assembly/components.bla', 'int'); for c=1:size(components,1) if ~isempty(components{c}) components{c}=components{c}+1; end; end;
if exist(sprintf('../aligns/%s_stage34.mat',curFNpref))==2
    eval(sprintf('load ../aligns/%s_stage34.mat %s_aligns_pa %s_pairs_asp; aligns_pa=%s_aligns_pa; aligns_asp=%s_pairs_asp; clear %s;',curFNpref,curFNpref,curFNpref,curFNpref,curFNpref,curFNpref));
else
    aligns_asp=[]; aligns_pa=[];
end

if size(specIndices,2)>1 specPositions = specIndices(:,2:3); specIndices=specIndices(:,1); else specPositions=zeros(length(specIndices),2); end;   numSpecs = size(specIndices,1);
fid = fopen(reportFN,'w'); if fid<=0 fprintf(1,'Error opening %s!\n',reportFN); return; end;

% Info about aligns
if isempty(aligns_asp) counts_asp=0; inSet_asp = 0; else
	matches_asp = zeros(size(aligns_asp,1),2);
	idx = csps_matchCols(aligns_asp(:,1),specIndices);   matches_asp(idx,1)=1;
	idx = csps_matchCols(aligns_asp(:,2),specIndices);   matches_asp(idx,2)=1;  
    inSet_asp = size(find(matches_asp(:,1)==1 & matches_asp(:,2)==1),1);
    counts_asp = size(find(matches_asp(:,1)==1 | matches_asp(:,2)==1),1)-inSet_asp;
end
if isempty(aligns_pa) counts_pa=0; inSet_pa=0; else
	matches_pa = zeros(size(aligns_pa,1),2);
	idx = csps_matchCols(aligns_pa(:,1),specIndices);   matches_pa(idx,1)=1;
	idx = csps_matchCols(aligns_pa(:,2),specIndices);   matches_pa(idx,2)=1;    
    inSet_pa = size(find(matches_pa(:,1)==1 & matches_pa(:,2)==1),1);
    counts_pa = size(find(matches_pa(:,1)==1 | matches_pa(:,2)==1),1)-inSet_pa;
end
fprintf(fid,'Number of spectra in set: %.0f\nProtein index: %.0f\n\nAlignments information:\n',size(specIndices,1),curProtIdx);
fprintf(fid,'  Matches within the set : %.0f (asp), %.0f (pa)\n',inSet_asp, inSet_pa);
fprintf(fid,'  Matches outside the set: %.0f (asp), %.0f (pa)\n\n',counts_asp, counts_pa);

% Info about components. Component types: matched to curProtIdx (1), matched to other prot (2), no reasonable tag (3), no spectrum from specIndices but contig seq still matches curProtIdx (4), no match to specIndices nor curProtIdx (0)
szComps = size(components,1);   compsInfo = zeros(szComps,5); % Component type (col.1), Number of spectra from specIndices (col.2), Total # spectra in comp (col.3), Index in tags (col.4), Protein position of the leftmost assembled spectrum (col.5)
szTags = size(tags,1);   overallMatched = [];   specMembership = zeros(max(specIndices),1); % Keep track of which component a spectrum belongs to
for compIdx=1:size(components,1)
    compsInfo(compIdx,3) = size(components{compIdx},2);   tagIdx = min(find([tags{:,1}]==compIdx));
    [common, idxSpecIndices] = intersect(specIndices',components{compIdx});   
    if isempty(common) 
        if ~isempty(tagIdx) & tagsMatch(tagIdx)==curProtIdx compsInfo(compIdx,1)=4; end; continue; 
    else
        overallMatched=[overallMatched common];
        specMembership(common) = compIdx;
        compsInfo(compIdx,2) = size(common,2);
        compsInfo(compIdx,5) = min(specPositions(idxSpecIndices,1));
        if compsInfo(compIdx,5)==0 & ~isempty(tagIdx) & size(tagsMatch,2)>1 compsInfo(compIdx,5)=tagsMatch(tagIdx,2); end;
    end;   
%     while curTag<=szTags & tags{curTag,1}<compIdx curTag=curTag+1; end;
    if isempty(tagIdx) compsInfo(compIdx,1)=3; continue; end;  % No reasonable tag found (according to tags)
    if tagsMatch(tagIdx)==curProtIdx compsInfo(compIdx,1)=1; else compsInfo(compIdx,1)=2; end;
    compsInfo(compIdx,4)=tagIdx;
end
specMembership = specMembership(specIndices);
comps_indices = {find(compsInfo(:,1)==1); find(compsInfo(:,1)==2); find(compsInfo(:,1)==3); find(compsInfo(:,1)==4);};
fprintf(fid,'Number of components with at least one spectrum from the set or matching the protein: %.0f\n',size(find(compsInfo(:,2)>0),1));
fprintf(fid,'   Tag matches same protein  : %.0f\n',size(comps_indices{1},1));
fprintf(fid,'   Tag matches other protein : %.0f\n',size(comps_indices{2},1));
fprintf(fid,'   No reasonable tag         : %.0f\n',size(comps_indices{3},1));
fprintf(fid,'   Matches protein, not specs: %.0f\n',size(comps_indices{4},1));
fprintf(fid,'Total # spectra in some component: %.0f (%.2f%%), %.0f not grouped\n\n',size(overallMatched,2),100*size(overallMatched,2)/size(specIndices,1),size(specIndices,1)-size(overallMatched,2));

% Info about identifications
comps_titles = {'*** Tag matches same protein ***';'*** Tag matches other protein ***';'*** No reasonable tag ***';'*** Tag matches protein, assembled spectra do not ***'};
for compType=1:4
    fprintf(fid,'Component type: %s\n----------------------------------------------------\n\n',comps_titles{compType});
    theseComps = [comps_indices{compType} compsInfo(comps_indices{compType},:)];
    if isempty(theseComps) continue; end;
    [foo,idxS] = sort(theseComps(:,6));   theseComps = theseComps(idxS,:);
    for cur=1:size(theseComps,1)
        compIdx = theseComps(cur,1);
        fprintf(fid,'Component %.0f, tags index %.0f, match starting at aa %.0f\n',compIdx,theseComps(cur,5),theseComps(cur,6));
        fprintf(fid,'   Spectra from set: %.0f (%.2f%%), %.0f others\n',theseComps(cur,3),100*theseComps(cur,3)/theseComps(cur,4),theseComps(cur,4)-theseComps(cur,3));
        fprintf(fid,'   DB match (0 mods): %s (ei=%.2f, ep=%.2f, tp=%.0f, fp=%.0f)\n',matchmaReport{compIdx}{1,8},matchmaReport{compIdx}{1,9},matchmaReport{compIdx}{1,10},matchmaReport{compIdx}{1,11},matchmaReport{compIdx}{1,12});
        if size(matchmaReport{compIdx},1)>1 fprintf(fid,'   DB match (1 mods): %s (ei=%.2f, ep=%.2f, tp=%.0f, fp=%.0f)\n',matchmaReport{compIdx}{2,8},matchmaReport{compIdx}{2,9},matchmaReport{compIdx}{2,10},matchmaReport{compIdx}{2,11},matchmaReport{compIdx}{2,12}); end;
        fprintf(fid,'   De-novo sequence : %s\n\n',matchmaReport{compIdx}{1,13});
    end;
    fprintf(fid,'\n');
end
fprintf(fid,'\n');

% Info about spectra
specStats = zeros(numSpecs,6); % sum(b,y expint) (ctgs/cPAe cols.1,2), bypercs: ctgs (cols 3,4) cPAe (cols.5,6)
if ~isempty(ctgs_expint)  specStats(:,1)=sum(ctgs_expint(specIndices,1:2)')'; end;
if ~isempty(ctgs_bypercs) specStats(:,3:4)=ctgs_bypercs(specIndices,[1 3]); end;
if ~isempty(cPAe_expint)  specStats(:,2)=sum(cPAe_expint(specIndices,1:2)')'; end;
if ~isempty(cPAe_bypercs) specStats(:,5:6)=cPAe_bypercs(specIndices,[1 3]); end;
specMatched = zeros(numSpecs,1);   specMatched(csps_matchCols(specIndices,overallMatched'))=1;
peptides = allPeptides(specIndices);
if ~isempty(specPositions)
    posOffset = min(min(specPositions));
    [specPositions, idxS] = sortrows(specPositions-posOffset+1);  % Left-align all positioned peptides
    specStats = specStats(idxS,:);   specMatched = specMatched(idxS,:);   peptides = peptides(idxS,:);   
    specMembership = specMembership(idxS,:);    specIndices = specIndices(idxS,:);
end; 
fprintf(fid,'Information about spectra in the set:\n');
for specIdx=1:numSpecs
    if specMatched(specIdx)==1 matchedStr=sprintf('%3.0f',specMembership(specIdx)); else matchedStr='   '; end;
    fprintf(fid,'%s (%.2f/%.2f) (b=%.2f/%.2f, y=%.2f/%.2f):',matchedStr,specStats(specIdx,1),specStats(specIdx,2),specStats(specIdx,3),specStats(specIdx,5),specStats(specIdx,4),specStats(specIdx,6));
    if isempty(specPositions) fprintf(fid,' %s (idx=%.0f)',peptides{specIdx},specIndices(specIdx)); else fprintf(fid,'%s[%4.0f] %s (idx=%.0f)',repmat(' ',1,specPositions(specIdx,1)),posOffset-1+specPositions(specIdx,1),peptides{specIdx},specIndices(specIdx)); end;
    if ~isempty(filenames) fprintf(fid, ' %s\n',filenames{specIndices(specIdx)}); else fprintf(fid, '\n'); end;
end
fclose(fid);
