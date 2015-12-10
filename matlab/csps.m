% Copyright 2007, The Regents of the University of California
% All Rights Reserved
%
% Permission to use, copy, modify and distribute any part of this
% program for educational, research and non-profit purposes, without fee,
% and without a written agreement is hereby granted, provided that the
% above copyright notice, this paragraph and the following three paragraphs
% appear in all copies.
%
% Those desiring to incorporate this work into commercial
% products or use for commercial purposes should contact the Technology
% Transfer & Intellectual Property Services, University of California,
% San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910,
% Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.
%
% IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
% FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
% INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN
% IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY
% OF SUCH DAMAGE.
%
% THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY
% OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
% ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO
% REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
% EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
% THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.

%
%  csps.m - top-level script to run Comparative Shotgun Protein Sequencing analysis
%
function res = csps(paramsFN, stage_str_cl)
% function res = csps(spsInfo,peakTol,pmTol)
%

if nargin<1
    fprintf(1,'Syntax: csps parameters_file_name\n\n');
    res=0; return;
end
if nargin<2 stage_str_cl=''; end

%
% Reading parameters
%
[status, params] = sn_read_params_file(paramsFN,1);  if status==0 res=-1; return; end;
mandatory_params = {'EXE_DIR';'CLUSTALW_EXE_DIR';'SPS_PROJECTS'};    params_names = fieldnames(params);
for idx=1:size(mandatory_params,1)
    if max(strcmp(mandatory_params{idx},params_names))==0 | isempty(getfield(params,mandatory_params{idx}))
        fprintf(1,'ERROR: Parameter %s is undefined in %s\n',mandatory_params{idx},paramsFN);    res=-1; return;
    end;
end;
pmTol = str2num(params.TOLERANCE_PM);    peakTol = str2num(params.TOLERANCE_PEAK);

tagLen    = str2num(params.TAG_LEN); % = 8;
numAAJmps = str2num(params.DOUBLE_AA_JUMPS); % = 1;
maxAAJump = str2num(params.MAX_AA_JUMP); % = 2;
tagFlanks = str2num(params.MATCH_TAG_FLANKING_MASSES); % = 0;

maxModMass = str2num(params.MAX_MOD_MASS);  % 100
minModMass = str2num(params.MIN_MOD_MASS);  % -100

v=version;   v = str2num(v(1));  if v>6 v='-v6'; else v=''; end;
STAGE_TAGSEARCH = 1;
STAGE_PROTID    = 2;
STAGE_MATCHMA   = 3;
STAGE_CLUSTALW  = 4;
STAGE_HOMGLUE   = 5;
STAGE_REPORT    = 6;
stage_str = lower(params.INITIAL_STAGE);    initial_stage=0;
if ~isempty(stage_str_cl) stage_str = stage_str_cl; end;
if strcmp(stage_str,'tagsearch')       initial_stage=STAGE_TAGSEARCH; end;
if strcmp(stage_str,'protid')          initial_stage=STAGE_PROTID; end;
if strcmp(stage_str,'alignment')       initial_stage=STAGE_MATCHMA; end;
if strcmp(stage_str,'clustalw')        initial_stage=STAGE_CLUSTALW; end;
if strcmp(stage_str,'sequencing')      initial_stage=STAGE_HOMGLUE; end;
if strcmp(stage_str,'report')          initial_stage=STAGE_REPORT; end;

%
% Initialize environment
%
c=clock;
curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
fprintf(1,'Comparative Shotgun Protein Sequencing 0.1.0: session started on %s %s\n',curDate,curTime);
dirContents = dir('.'); dirContents = {dirContents(:).name}';
if isempty(find(strcmp('spectra',dirContents))) status = dos('mkdir spectra'); end;
if isempty(find(strcmp('clustalw',dirContents))) status = dos('mkdir clustalw'); end;
if isempty(find(strcmp('homology',dirContents))) status = dos('mkdir homology'); end;
if isempty(find(strcmp('report',dirContents))) status = dos('mkdir report'); end;

global AAletters AAmasses AAcodes AAnames;
global bProfile yProfile bNames yNames;
load(sprintf('%s/msparams.mat',params.EXE_DIR));
if ~isempty(params.AMINO_ACID_MASSES) [res,AAmasses,AAletters] = sn_read_masses(params.AMINO_ACID_MASSES,AAmasses,AAletters); if res<=0 fprintf(1,'Error reading amino acid masses from: %s\n',params.AMINO_ACID_MASSES); res=-1; return; end; end;


% Read SPS projects information
spsLines = sn_load_lines(params.SPS_PROJECTS);   if isempty(spsLines) fprintf(1,'ERROR: Could not load %s!\n',params.SPS_PROJECTS); res=-1; return; end;
numProjs = size(spsLines,1);   spsInfo=cell(numProjs,4);   pIdx=1;
for lineIdx=1:numProjs
    if ~isempty(spsLines{lineIdx}) & spsLines{lineIdx}(1)~='#'
        idxSep = find(spsLines{lineIdx}==';');
        if isempty(idxSep) | max(idxSep)~=length(spsLines{lineIdx}) idxSep=[0 idxSep length(spsLines{lineIdx})+1]; end;
        if length(idxSep)~=5 fprintf(1,'ERROR: Invalid line %d in %s (%s)!\n',lineIdx,params.SPS_PROJECTS,spsLines{lineIdx}); res=-1; return; end;

        spsInfo{pIdx,1} = spsLines{lineIdx}( (idxSep(1)+1):(idxSep(2)-1) );
        spsInfo{pIdx,2} = spsLines{lineIdx}( (idxSep(2)+1):(idxSep(3)-1) );
        spsInfo{pIdx,3} = str2num(spsLines{lineIdx}( (idxSep(3)+1):(idxSep(4)-1) ));
        spsInfo{pIdx,4} = str2num(spsLines{lineIdx}( (idxSep(4)+1):(idxSep(5)-1) ));
        pIdx=pIdx+1;
    end;
end
spsInfo=spsInfo(1:pIdx-1,:);   if isempty(spsInfo) fprintf(1,'ERROR: Could not load %s!\n',params.SPS_PROJECTS); res=-1; return; end;

%
% Tag generation + database search
%
cd('spectra');
if initial_stage <= STAGE_TAGSEARCH
	c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
	curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
    fprintf(1,'Starting stage: tagsearch (%s %s)\n',curDate,curTime);

    stgParams = {};
    if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../%s',params.AMINO_ACID_MASSES)}]; end;
    for projIdx=1:numProjs
        fid = fopen(sprintf('tagsearch_%s.params',spsInfo{projIdx,1}),'w');   if fid<=0 fprintf(1,'ERROR opening tagsearch_%s.params!\n',spsInfo{projIdx,1}); res=-1; return; end;
		fprintf(fid,'%d\nINPUT_SPECS_PKLBIN=%s/assembly/sps_seqs.pklbin\nINPUT_FASTA=../%s\n',10+size(stgParams,1),spsInfo{projIdx,2},params.FASTA_DATABASE);
		fprintf(fid,'TOLERANCE_PEAK=%f\nTOLERANCE_PM=%f\nTAG_LEN=%d\nDOUBLE_AA_JUMPS=%d\nMATCH_TAG_FLANKING_MASSES=%d\nMAX_PARSIMONY=0\n',spsInfo{projIdx,3},spsInfo{projIdx,4},tagLen,numAAJmps,tagFlanks);
		fprintf(fid,'OUTPUT_PEPTIDES=tagsearch_%s.txt\nOUTPUT_MATCHED_PROTS=tagsearch_%s.bla\n',spsInfo{projIdx,1},spsInfo{projIdx,1});
		for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
		fclose(fid);
		status = dos(sprintf('%s/tagsearch tagsearch_%s.params',params.EXE_DIR,spsInfo{projIdx,1}));    if status~=0 fprintf(1,'ERROR executing %s/tagsearch!\n',params.EXE_DIR); res=-1; return; end;
    end
end;
cd('..');


% Run protid
if initial_stage<=STAGE_PROTID
	c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
	curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
    fprintf(1,'Starting stage: protid (%s %s)\n',curDate,curTime);

    stgParams = {};
    if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=%s',params.AMINO_ACID_MASSES)}]; end;

    fid = fopen('protid_list.txt','w');   if fid<=0 fprintf(1,'ERROR opening protid_list.txt!\n'); res=-1; return; end;
    for projIdx=1:numProjs fprintf(fid,'spectra/tagsearch_%s.bla\n',spsInfo{projIdx,1}); end
	fclose(fid);
    fid = fopen('protid.params','w');   if fid<=0 fprintf(1,'ERROR opening protid.params!\n'); res=-1; return; end;
    fprintf(fid,'%d\nINPUT_FASTA=%s\nOUTPUT_FASTA=protid.fasta\nINPUT_MATCHES_LIST=protid_list.txt\n',3+size(stgParams,1),params.FASTA_DATABASE);
    for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
	fclose(fid);
	status = dos(sprintf('%s/protid protid.params',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/protid!\n',params.EXE_DIR); res=-1; return; end;
    delete('protid_list.txt');
end
protid_db = csps_load_fasta('protid.fasta','');

% Run matchma against reduced database
cd('homology');
if initial_stage<=STAGE_MATCHMA
    mmaData = cell(numProjs,6);  % matchma output (i,:): #matches (1), matched protein info (2), selected b/y spectra (3), matched indices (4), indices of kept contigs (5), report (6)
    stgParams = {sprintf('MIN_NUM_MATCH_PEAKS=%s',params.MIN_MATCHED_PEAKS_DB)};
    if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../%s',params.AMINO_ACID_MASSES)}]; end;

    for projIdx=1:numProjs
        fid=fopen(sprintf('matchma_%s.params',spsInfo{projIdx,1}),'w'); if fid<=0 fprintf(1,'ERROR opening matchma_%s.params!\n',spsInfo{projIdx,1}); res=-1; return; end;
        fprintf(fid,'%d\nINPUT_CLUSTERS=%s/assembly/path_spectra_as_cluster.txt\nINPUT_FASTA=../protid.fasta\nOUTPUT_CSV=matchma_%s_csv.txt\n',15+size(stgParams,1),spsInfo{projIdx,2},spsInfo{projIdx,1});
        fprintf(fid,'TOLERANCE_PEAK=%f\nTOLERANCE_PM=%f\nENFORCE_ENDPEAKS=0\n',spsInfo{projIdx,3},spsInfo{projIdx,4});
        fprintf(fid,'RESOLUTION=0.1\nMIN_RATIO=0.4\nSPEC_TYPE_MSMS=0\nMIN_MOD_MASS=%f\nMAX_MOD_MASS=%f\nMAX_NUM_MODS=%s\n',minModMass,maxModMass,params.MAX_NUM_MODS);
        fprintf(fid,'OUTPUT_MATCHED_PROTS=matchma_%s_mp.bin\nOUTPUT_MATCHED_SPECS=matchma_%s_mspecs.pklbin\nOUTPUT_MATCHED_PEAKS_IDX=matchma_%s_midx.pklbin\n',spsInfo{projIdx,1},spsInfo{projIdx,1},spsInfo{projIdx,1});
        for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
        fclose(fid);
        status = dos(sprintf('%s/matchma matchma_%s.params',params.EXE_DIR,spsInfo{projIdx,1}));    if status~=0 fprintf(1,'ERROR executing %s/matchma matchma_%s.params!\n',params.EXE_DIR,spsInfo{projIdx,1}); res=-1; return; end;

        mmaData{projIdx,2} = sn_load_binarray(sprintf('matchma_%s_mp.bin',spsInfo{projIdx,1}),'int');
        mmaData{projIdx,3} = sn_load_pklbin(sprintf('matchma_%s_mspecs.pklbin',spsInfo{projIdx,1}));
        mmaData{projIdx,4} = sn_load_pklbin(sprintf('matchma_%s_midx.pklbin',spsInfo{projIdx,1}));
        
        % Select only good matches
%         masab_specs = sn_load_pklbin(sprintf('%s/assembly/sps_seqs.pklbin',spsInfo{projIdx,2}));
%         [report,reportScores]=csps_report_matchma(masab_specs,sprintf('matchma_%s_csv.txt',spsInfo{projIdx,1}), [], maxModMass, peakTol,str2num(params.MAX_NUM_MODS)+1,sprintf('matchma_%s_sps_report.html',spsInfo{projIdx,1}));
%         tags = csps_filterLongestTag(report, str2num(params.REPORT_MIN_TAG), maxAAJump, maxModMass);
%         fid=fopen(sprintf('matchma_%s_tags.txt',spsInfo{projIdx,1}),'w'); if fid<=0 fprintf(1,'ERROR opening matchma_%s_tags.txt!\n',spsInfo{projIdx,1}); res=-1; return; end;
%         for tagIdx=1:size(tags,1)
%             fprintf(fid,'%d;[%.1f]%s[%.1f]\n',tags{tagIdx,1},tags{tagIdx,4},tags{tagIdx,2},tags{tagIdx,5});
%         end
%         fclose(fid);
%         if ~isempty(v) save('tags','-v6','tags'); else save('tags','tags'); end;
% 
%         mmaData{projIdx,5} = [tags{:,1}]';
        mmaData{projIdx,5} = find(mmaData{projIdx,2}(:,1)>=0);
        mmaData{projIdx,1} = length(mmaData{projIdx,5});
        mmaData{projIdx,2} = mmaData{projIdx,2}(mmaData{projIdx,5},:);
        mmaData{projIdx,3} = mmaData{projIdx,3}(mmaData{projIdx,5},:);
        mmaData{projIdx,4} = mmaData{projIdx,4}(mmaData{projIdx,5},:);
    end
    
    % Merge matchma results
    numContigs = sum([mmaData{:,1}]);   mma_protInfo = zeros(numContigs,3);   mma_specs = cell(numContigs,5);   mma_mIdx = cell(numContigs,5);   
    mma_index = zeros(numContigs,2);    % Index of retained contigs: project index (col.1), contig index (col.2)
    idx = [0 cumsum([mmaData{:,1}])];
    for projIdx=1:numProjs
        range = [idx(projIdx)+1:idx(projIdx+1)];    mma_protInfo(range,:) = mmaData{projIdx,2};
        mma_specs(range,:) = mmaData{projIdx,3};    mma_mIdx(range,:) = mmaData{projIdx,4};
        mma_index(range,:) = [projIdx*ones(mmaData{projIdx,1},1) mmaData{projIdx,5}];
    end
    sn_save_pklbin(mma_specs,'../spectra/contigs.pklbin');
    sn_save_pklbin(mma_mIdx,'../spectra/contigs_matches.pklbin');
    sn_save_binarray('../spectra/contigs_matches.bin',mma_protInfo,'uint');
    sn_save_binarray('../spectra/contigs_index.bin',mma_index,'uint');
    if ~isempty(v) save('mmaData','-v6','mmaData'); else save('mmaData','mmaData'); end;
else
    mma_specs = sn_load_pklbin('../spectra/contigs.pklbin');
    mma_mIdx  = sn_load_pklbin('../spectra/contigs_matches.pklbin');
    mma_protInfo = sn_load_binarray('../spectra/contigs_matches.bin','uint');
    mma_index    = sn_load_binarray('../spectra/contigs_index.bin','uint');
    load mmaData;
end

% Running clustalw testLC.fasta -> testLC.aln outputs text below ; keep only alignments with score >= 250 (parameterize)
cd('../clustalw');
% Determine master protein - that with the largest number of matched masses
% matchma output (i,:): #matches (1), matched protein info (2), selected b/y spectra (3), matched indices (4),
%                          indices of kept contigs (5), report (6)
dbSize = size(protid_db,1);   numMatches = zeros(dbSize,1);
for projIdx=1:numProjs
    for matchIdx=1:mmaData{projIdx,1}
        pIdx = mmaData{projIdx,2}(matchIdx,1)+1;
fprintf(1,'pIdx = %d, dbSize = %d, size(numMatches,1) = %d\n',pIdx,dbSize, size(numMatches,1));
fprintf(1,'size(mmaData) = (%d,%d), projIdx = %d, matchIdx = %d\n',size(mmaData,1),size(mmaData,2),projIdx,matchIdx);
fprintf(1,'size(mmaData{projIdx,4},1) = %d\n',size(mmaData{projIdx,4},1));
fprintf(1,'size(mmaData{projIdx,4}{matchIdx,2},1) = %d\n',size(mmaData{projIdx,4}{matchIdx,2},1));
        numMatches(pIdx) = numMatches(pIdx)+size(mmaData{projIdx,4}{matchIdx,2},1);
    end
end
[foo,mIdx] = max(numMatches);   if length(mIdx)>1 mIdx=mIdx(1); end;
if str2num(params.FORCE_REFERENCE)>0 
    mIdx=str2num(params.FORCE_REFERENCE);
    fprintf(1,'(user forced reference): ');
end
fprintf(1,'Main matched protein %d: %s %s\n',mIdx,protid_db{mIdx,1},protid_db{mIdx,3});

if initial_stage<=STAGE_CLUSTALW & ~isempty(params.CLUSTALW_EXE_DIR)
    % Compute clustalw alignments between master protein and all other proteins
    numHomologProts = 0;
    fidCWIndex = fopen(sprintf('cwindex.txt',mIdx,pIdx),'w');  if fidCWIndex<=0 fprintf(1,'ERROR creating cwindex.txt!\n'); res=-1; return; end;
    for pIdx=1:dbSize
        if pIdx~=mIdx
            fid=fopen(sprintf('cwseqs_%d_%d.fasta',mIdx,pIdx),'w');  if fid<=0 fprintf(1,'ERROR creating cwseqs_%d_%d.fasta!\n',mIdx,pIdx); res=-1; return; end;
            fprintf(fid,'>Protein_%d\n%s\n\n>Protein_%d\n%s\n\n',mIdx,protid_db{mIdx,2},pIdx,protid_db{pIdx,2}); 
            fclose(fid);
            [foo,txt] = dos(sprintf('%s/clustalw cwseqs_%d_%d.fasta',params.CLUSTALW_EXE_DIR,mIdx,pIdx));
            txt=txt(max(findstr(txt,'Alignment Score '))+16:length(txt));    txt=txt(1:min(find(txt==10))-1);
            if str2num(txt)>=str2num(params.CLUSTALW_MINSCORE)
                fprintf(fidCWIndex,'%d;%d;../clustalw/cwseqs_%d_%d.aln\n',mIdx,pIdx,mIdx,pIdx); % Output zero-based protein indices
                numHomologProts = numHomologProts+1;
            end
%             cwScores(pIdx) = str2num(txt);
%             [cwAligns{pIdx,1}, cwAligns{pIdx,2}, cwAligns{pIdx,3}] = clustalw_load_aln(sprintf('cwseqs_%d_%d.fasta',mIdx,pIdx));
        end
    end
    fclose(fidCWIndex);   if ~isempty(v) save('cwdata','-v6','numHomologProts'); else save('cwdata','numHomologProts'); end;
else load cwdata; end;
cd('../homology');

% -----------------------------------------------------
% 
% 
% 
%  CLUSTAL W (1.83) Multiple Sequence Alignments
% 
% 
% 
% Sequence format is Pearson
% Sequence 1: Sequence                    219 aa
% Sequence 2: gi|42543442|pdb|1Q9V|A      219 aa
% Start of Pairwise alignments
% Aligning...
% Sequences (1:2) Aligned. Score:  97
% Guide tree        file created:   [testLC.dnd]
% Start of Multiple Alignment
% There are 1 groups
% Aligning...
% Group 1: Sequences:   2      Score:4702
% Alignment Score 1291
% CLUSTAL-Alignment file created  [testLC.aln]
% -----------------------------------------------------

% Run homglue
spsNames = cell(size(mma_index,1),1);   for i=1:size(spsNames,1) spsNames{i}=sprintf('%s:%d',spsInfo{mma_index(i,1),1},mma_index(i,2)); end;
sn_save_lines('ref_sps_names.txt',spsNames,1);
if initial_stage<=STAGE_HOMGLUE
    stgParams={};
    if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../%s',params.AMINO_ACID_MASSES)}]; end;

    if numHomologProts>0 stgParams = [stgParams; {sprintf('INPUT_HOMOLOGIES=../clustalw/cwindex.txt\n')}]; end;
	stgParams = [stgParams; {sprintf('MAX_MOD_MASS=%s',params.MAX_MOD_MASS)}];
    fid=fopen('homglue.params','w'); if fid<=0 fprintf(1,'ERROR opening homglue.params!\n'); res=-1; return; end;
    fprintf(fid,'%d\nINPUT_SPECS_PKLBIN=../spectra/contigs.pklbin\nINPUT_FASTA=../protid.fasta\n',14+size(stgParams,1));
    fprintf(fid,'INPUT_MATCHED_PEAKS_IDX=../spectra/contigs_matches.pklbin\nINPUT_MATCHED_PROTS=../spectra/contigs_matches.bin\n');
    fprintf(fid,'OUTPUT_CSV=homglue_matches.txt\nOUTPUT_SPECS=homglue_matches.pklbin\nOUTPUT_MATCHES_REF=homglue_ref\nOUTPUT_MATCHES_CSPS=homglue_matches\n');
    fprintf(fid,'TOLERANCE_PEAK=%f\nTOLERANCE_PM=%f\n',peakTol,pmTol);
    fprintf(fid,'RESOLUTION=0.1\nSPEC_TYPE_MSMS=0\nGRAPH_TYPE=2\nMIN_CONTIG_SET=1\n');
    for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
    fclose(fid);
    status = dos(sprintf('%s/homglue homglue.params',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/homglue homglue.params!\n',params.EXE_DIR); res=-1; return; end;
    
    % No-glues reference alignments (for comparison purposes only)
    fid=fopen('homglue_ref.params','w'); if fid<=0 fprintf(1,'ERROR opening homglue_ref.params!\n'); res=-1; return; end;
    fprintf(fid,'%d\nINPUT_SPECS_PKLBIN=../spectra/contigs.pklbin\nINPUT_FASTA=../protid.fasta\n',14+size(stgParams,1));
    fprintf(fid,'INPUT_MATCHED_PEAKS_IDX=../spectra/contigs_matches.pklbin\nINPUT_MATCHED_PROTS=../spectra/contigs_matches.bin\n');
    fprintf(fid,'OUTPUT_CSV=homglue_ref.txt\nOUTPUT_SPECS=homglue_ref.pklbin\nSKIP_GLUES=1\n');
    fprintf(fid,'TOLERANCE_PEAK=%f\nTOLERANCE_PM=%f\nINPUT_SPECS_NAMES=ref_sps_names.txt\n',peakTol,pmTol);
    fprintf(fid,'RESOLUTION=0.1\nSPEC_TYPE_MSMS=0\nGRAPH_TYPE=2\nMIN_CONTIG_SET=1\n');
    for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
    fclose(fid);
    status = dos(sprintf('%s/homglue homglue_ref.params',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/homglue homglue.params!\n',params.EXE_DIR); res=-1; return; end;
end
% cd('..');

% Report homglue results
homglue_specs = sn_load_pklbin('homglue_matches.pklbin');
[matchStats, protStats, matchParcels] = csps_report('main','homglue_matches.txt',size(homglue_specs,1),protid_db,mIdx,0,'','',peakTol,20);
if ~isempty(v) save('homglue_results',v,'matchStats','protStats','matchParcels'); else save('homglue_results','matchStats','protStats','matchParcels'); end;

homglue_ref_specs = sn_load_pklbin('homglue_ref.pklbin');
[matchStats, protStats, matchParcels] = csps_report('ref','homglue_ref.txt',size(homglue_ref_specs,1),protid_db,mIdx,0,'','',peakTol,20,spsNames);
cd('..');

res=0;
