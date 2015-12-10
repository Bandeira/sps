% Copyright 2007,2008,2009 The Regents of the University of California
% All Rights Reserved
%
% Permission to use, copy, modify and distribute any part of this
% program for educational, research and non-profit purposes, by non-profit 
% institutions only, without fee, and without a written agreement is hereby 
% granted, provided that the above copyright notice, this paragraph and 
% the following three paragraphs appear in all copies.
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
%  specnets.m - top-level script to run Spectral Networks analysis
%
function res = specnets(paramsFN, stage_str_cl, runMode)

if nargin<1
    fprintf(1,'Syntax: specnets parameters_file_name\n\n');
    res=0; return;
end
if nargin<2 stage_str_cl=''; end
if nargin<3 runMode='specnets'; end

%
% Reading parameters
%
[status, params] = sn_read_params_file(paramsFN,1);  if status==0 res=-1; return; end;

mandatory_params = {'INPUT_SPECS_MS';'EXE_DIR'};    params_names = fieldnames(params);
for idx=1:size(mandatory_params,1)
    if max(strcmp(mandatory_params{idx},params_names))==0 | isempty(getfield(params,mandatory_params{idx}))
        fprintf(1,'ERROR: Parameter %s is undefined in %s\n',mandatory_params{idx},paramsFN);    res=-1; return;
    end;
end;
pmTol = str2num(params.TOLERANCE_PM);    peakTol = str2num(params.TOLERANCE_PEAK);   projDir = pwd;
if strcmp(runMode,'sps') pairsNumCols=6; else pairsNumCols=5; end;

STAGE_SCORING = 1;
STAGE_FILTERPAIRS = 2;
STAGE_ALIGNMENT = 3;
STAGE_FILTERSTARPAIRS = 4;
STAGE_ASSEMBLY = 5;
STAGE_TAGSEARCH = 6;
STAGE_SPECNETS = 7;
STAGE_REPORT = 8;
stage_str = lower(params.INITIAL_STAGE);    initial_stage=0;
if ~isempty(stage_str_cl) stage_str = stage_str_cl; end;
if strcmp(stage_str,'scoring')         initial_stage=STAGE_SCORING; end;
if strcmp(stage_str,'filterpairs')     initial_stage=STAGE_FILTERPAIRS; end;
if strcmp(stage_str,'alignment')       initial_stage=STAGE_ALIGNMENT; end;
if strcmp(stage_str,'filterstarpairs') initial_stage=STAGE_FILTERSTARPAIRS; end;
if strcmp(stage_str,'assembly')        initial_stage=STAGE_ASSEMBLY; end;
if strcmp(stage_str,'tagsearch')       initial_stage=STAGE_TAGSEARCH; end;
if strcmp(stage_str,'specnets')        initial_stage=STAGE_SPECNETS; end;
if strcmp(stage_str,'report')          initial_stage=STAGE_REPORT; end;

%
% Initialize environment
%
c=clock;
curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
fprintf(1,'Spectral Networks 2.0.0: session started on %s %s\n',curDate,curTime);
dirContents = dir('.'); dirContents = {dirContents(:).name}';
if isempty(find(strcmp('spectra',dirContents))) status = dos('mkdir spectra'); end;
if isempty(find(strcmp('aligns',dirContents))) status = dos('mkdir aligns'); end;
if isempty(find(strcmp('specnets',dirContents))) status = dos('mkdir specnets'); end;
if isempty(find(strcmp('assembly',dirContents))) status = dos('mkdir assembly'); end;
if isempty(find(strcmp('report',dirContents))) status = dos('mkdir report'); end;


global AAletters AAmasses AAcodes AAnames;
global bProfile yProfile bNames yNames;
load(sprintf('%s/msparams.mat',params.EXE_DIR));
if ~isempty(params.AMINO_ACID_MASSES) [res,AAmasses,AAletters] = sn_read_masses(params.AMINO_ACID_MASSES,AAmasses,AAletters); if res<=0 fprintf(1,'Error reading amino acid masses from: %s\n',params.AMINO_ACID_MASSES); res=-1; return; end; end;

%
% Loading MS/MS data (mgf/ms2/pkl with multiple spectra per file)
%
if initial_stage==0
    % Find out information about the files to load
    idxSep = find(params.INPUT_SPECS_MS==';');   if isempty(idxSep) | max(idxSep)~=length(params.INPUT_SPECS_MS) idxSep=[idxSep length(params.INPUT_SPECS_MS)+1]; end;
    numFiles = length(idxSep);   fileinfo = cell(numFiles,4); % complete filename, no extension (1), filetype (2), # MS/MS specs in file (3), project file name: specs_ms_%d or clusters_0_%d (4)
    idxSep = [0 idxSep];
    fid_index = fopen('spectra/input_index.txt','w'); if fid_index<=0 fprintf(1,'ERROR: Could not create input_index.txt\n'); res=-1; return; end;
    numSpecs = 0;
    for fIdx=1:numFiles
        fileinfo{fIdx,1} = params.INPUT_SPECS_MS( (idxSep(fIdx)+1):(idxSep(fIdx+1)-1) );

        % Find file type
        idx = max(find(fileinfo{fIdx,1}=='.'));
		if ~isempty(idx) fileinfo{fIdx,2}=fileinfo{fIdx,1}(idx+1:length(fileinfo{fIdx,1})); else fileinfo{fIdx,2}=''; end;
		if isempty(fileinfo{fIdx,2}) | (strcmp(lower(fileinfo{fIdx,2}),'mgf')==0 & strcmp(lower(fileinfo{fIdx,2}),'ms2')==0 & strcmp(lower(fileinfo{fIdx,2}),'mzxml')==0 & strcmp(lower(fileinfo{fIdx,2}),'pkl')==0 & strcmp(lower(fileinfo{fIdx,2}),'pklbin')==0)
            fprintf(1,'ERROR: Unknown file type %s (parameter INPUT_SPECS_MS in %s, file %s)\n',fileinfo{fIdx,2},paramsFN,fileinfo{fIdx,1}); res=-1; return;
		end;
        fileinfo{fIdx,1}=fileinfo{fIdx,1}(1:idx-1);

        % Convert to .pklbin and get number of spectra
        if(strcmp(lower(fileinfo{fIdx,2}),'pklbin'))
            specs_raw = sn_load_pklbin(sprintf('%s.%s',fileinfo{fIdx,1},fileinfo{fIdx,2}),1);  if isempty(specs_raw) fprintf(1,'ERROR reading %s.%s!\n',fileinfo{fIdx,1},fileinfo{fIdx,2}); res=-1; return; end;
            sn_save_pklbin(specs_raw, sprintf('spectra/specs_ms_%d.pklbin',fIdx));
        else
            status = -1;
            if(strcmp(lower(fileinfo{fIdx,2}),'mzxml'))  status = dos(sprintf('%s/convert mzxml %s.%s spectra/specs_ms_%d',params.EXE_DIR,fileinfo{fIdx,1},fileinfo{fIdx,2},fIdx)); end;
% status = 0;
            if(strcmp(lower(fileinfo{fIdx,2}),'ms2'))    status = dos(sprintf('%s/convert ms2   %s.%s spectra/specs_ms_%d',params.EXE_DIR,fileinfo{fIdx,1},fileinfo{fIdx,2},fIdx)); end;
            if(strcmp(lower(fileinfo{fIdx,2}),'mgf'))    status = dos(sprintf('%s/convert mgf   %s.%s spectra/specs_ms_%d',params.EXE_DIR,fileinfo{fIdx,1},fileinfo{fIdx,2},fIdx)); end;
            if(strcmp(lower(fileinfo{fIdx,2}),'pkl'))    status = dos(sprintf('%s/convert pkl   %s.%s spectra/specs_ms_%d',params.EXE_DIR,fileinfo{fIdx,1},fileinfo{fIdx,2},fIdx)); end;
            if status~=0 fprintf(1,'ERROR reading %s.%s!\n',fileinfo{fIdx,1},fileinfo{fIdx,2}); res=-1; return; end;
        end;
        fid=fopen(sprintf('spectra/specs_ms_%d.pklbin',fIdx),'r'); v=fread(fid,2,'uint'); if v(1)>0 fileinfo{fIdx,3}=v(1); else fileinfo{fIdx,3}=v(2); end; fclose(fid);
        numSpecs = numSpecs + fileinfo{fIdx,3};
        fileinfo{fIdx,4} = sprintf('specs_ms_%d.pklbin',fIdx);
        fprintf(1,' - Read %s.%s: %d spectra\n',fileinfo{fIdx,1},fileinfo{fIdx,2},fileinfo{fIdx,3});
        fprintf(fid_index,'%s.%s\n',fileinfo{fIdx,1},fileinfo{fIdx,2});
    end;
    fprintf(1,'Total: %d spectra\n',numSpecs);
    fileinfoMSMS = fileinfo;   fclose(fid_index);

    cd('spectra');
    if str2num(params.CLUSTER_MIN_SIZE)>0
        dirContents = dir('.'); dirContents = {dirContents(:).name}';
        if isempty(find(strcmp('tmp',dirContents))) status = dos('mkdir tmp'); end;
        if isempty(find(strcmp('out',dirContents))) status = dos('mkdir out'); end;

        fid=fopen('cluster_files.txt','w'); if fid<=0 fprintf(1,'ERROR opening cluster_files.txt!\n'); res=-1; return; end;
        for fIdx=1:numFiles
            cur_specs = sn_load_pklbin(sprintf('specs_ms_%d.pklbin',fIdx),1);  if isempty(cur_specs) fprintf(1,'ERROR reading specs_ms_%d.pklbin!\n',fIdx); res=-1; return; end;
% Preprocess spectra here (if needed: e.g. fake z>3 as z=3)
%             for s=1:size(cur_specs,1) 
%                 if ~isempty(cur_specs{s,4}) 
%                     if cur_specs{s,4}>3 cur_specs{s,4}=3; end;
%                 end;
%             end;
%            delete(sprintf('specs_ms_%d.pklbin',fIdx));
            sn_save_mgf(cur_specs(:,2), cur_specs(:,3), cur_specs(:,4), sprintf('specs_ms_%d.mgf',fIdx));
            fprintf(fid,'specs_ms_%d.mgf\n',fIdx);
        end
        fclose(fid);
        
%
%         Params with MS-Cluster syntax before March 15, 2010
%
%         if ~isempty(params.MIN_SPECTRUM_QUALITY) | strcmp(lower(params.GUESS_CHARGE),'yes') 
%             clusterMode = sprintf('-filter_model_name %s',params.CLUSTER_MODEL);
%             if ~isempty(params.MIN_SPECTRUM_QUALITY) clusterMode=strcat(clusterMode,sprintf(' -min_filter_prob %s',params.MIN_SPECTRUM_QUALITY)); end
%             if strcmp(lower(params.GUESS_CHARGE),'yes') clusterMode=strcat(clusterMode,' -assign_charges'); end;
%         else clusterMode=''; end
%         [status,output] = dos(sprintf('%s/MSCluster_bin -list cluster_files.txt -min_size %s -name clusters -slice_width %.3f -model_dir %s/Models_mscluster %s',params.EXE_DIR,params.CLUSTER_MIN_SIZE,2*pmTol,params.EXE_DIR,clusterMode),'-echo'); if status~=0 fprintf(1,'ERROR executing %s/MSCluster_bin!\n',params.EXE_DIR); res=-1; return; end;

%         if ~isempty(params.MIN_SPECTRUM_QUALITY) | strcmp(lower(params.GUESS_CHARGE),'yes') 
%             clusterMode = sprintf('--model %s',params.CLUSTER_MODEL);
            clusterMode = sprintf('--assign-charges --model %s',params.CLUSTER_MODEL);
            if ~isempty(params.MIN_SPECTRUM_QUALITY) clusterMode=strcat(clusterMode,sprintf(' --sqs %s',params.MIN_SPECTRUM_QUALITY)); end
%             if strcmp(lower(params.GUESS_CHARGE),'yes') clusterMode=strcat(clusterMode,' --assign-charges'); end;
%         else clusterMode=''; end
        if(~isempty(params.CLUSTER_PMTOL_PPM))
            clusterMode=strcat(clusterMode,sprintf(' --precursor-ppm %s',params.CLUSTER_PMTOL_PPM));
        end;
        msclusterCmd = sprintf('%s/MsCluster_bin --list cluster_files.txt --output-name clusters --window %.3f --model-dir %s/Models_mscluster --fragment-tolerance %s %s',params.EXE_DIR,2*pmTol,params.EXE_DIR,params.TOLERANCE_PEAK,clusterMode);
        fprintf(1,msclusterCmd);
        [status,output] = dos(msclusterCmd,'-echo'); if status~=0 fprintf(1,'ERROR executing %s/MSCluster_bin!\n',params.EXE_DIR); res=-1; return; end;
% output = '';
        sn_save_lines('../log_mscluster.txt',{output},0);

%         for fIdx=1:numFiles delete(sprintf('specs_ms_%d.mgf',fIdx)); end
        cd('tmp');   delete('*.*');   cd('..');   rmdir('tmp');

%         cd('out');   dirContents = dir('./*.mgf');    dirContents = {dirContents(:).name}';   numFiles = size(dirContents,1);
%         fileinfo = cell(numFiles,4);
%         for cIdx=1:numFiles
%             fileinfo{cIdx,1} = sprintf('clusters_0_%d',cIdx);   fileinfo{cIdx,2} = 'mgf';   fileinfo{cIdx,4} = sprintf('clusters_%d.pklbin',cIdx);
%             status = dos(sprintf('%s/convert mgf clusters_0_%d.mgf ../clusters_%d',params.EXE_DIR,cIdx,cIdx)); if status~=0 fprintf(1,'ERROR converting MSCluster output file clusters_0_%d.mgf!\n',cIdx); res=-1; return; end;
%             fid=fopen(sprintf('../clusters_%d.pklbin',cIdx),'r'); v=fread(fid,2,'uint'); if v(1)>0 fileinfo{cIdx,3}=v(1); else fileinfo{cIdx,3}=v(2); end; fclose(fid);
%         end;
%         cd('..');

        cd('out/mgf');   dirContents = dir('./*.mgf');    dirContents = {dirContents(:).name}';   numFiles = size(dirContents,1);
        fileinfo = cell(numFiles,4);
        for cIdx=1:numFiles
            fileinfo{cIdx,1} = sprintf('clusters_0_%d',cIdx);   fileinfo{cIdx,2} = 'mgf';   fileinfo{cIdx,4} = sprintf('clusters_%d.pklbin',cIdx);
            status = dos(sprintf('%s/convert mgf %s ../../clusters_%d',params.EXE_DIR,dirContents{cIdx},cIdx)); if status~=0 fprintf(1,'ERROR converting MSCluster output file %s!\n',dirContents{cIdx}); res=-1; return; end;
            fid=fopen(sprintf('../../clusters_%d.pklbin',cIdx),'r'); v=fread(fid,2,'uint'); if v(1)>0 fileinfo{cIdx,3}=v(1); else fileinfo{cIdx,3}=v(2); end; fclose(fid);
        end;
        cd('../..');
    end
    
    % Merge-load the input files
    numSpecs = [fileinfo{:,3}];   specs_raw = cell(sum(numSpecs),5);   filenames = cell(size(specs_raw,1),1);    numSpecs=[0 cumsum(numSpecs)];
    for fIdx=1:numFiles
        cur_specs = sn_load_pklbin(fileinfo{fIdx,4},1);  if isempty(cur_specs) fprintf(1,'ERROR reading %s!\n',fileinfo{fIdx,4}); res=-1; return; end;
        specs_raw(numSpecs(fIdx)+1:numSpecs(fIdx+1),:) = cur_specs;   clear cur_specs;
        if strncmp(fileinfo{fIdx,4},'clusters_',9) delete(fileinfo{fIdx,4}); end;
%         for specIdx=numSpecs(fIdx)+1:numSpecs(fIdx+1) filenames{specIdx}=sprintf('%s.%s:%d',fileinfo{fIdx,1},fileinfo{fIdx,2},specIdx-numSpecs(fIdx)); end;
    end
    delete('clusters_*.bin');
    for specIdx=1:size(filenames,1) filenames{specIdx}=sprintf('specs_ms.mgf:%d',specIdx); end;
    sn_save_pklbin(specs_raw, 'specs_ms.pklbin');   sn_save_lines('specs_ms.txt',filenames,1);
else
	cd('spectra');
	specs_raw = sn_load_pklbin('specs_ms.pklbin',1);  if isempty(specs_raw) fprintf(1,'ERROR reading specs_ms.pklbin!\n'); res=-1; return; end;
    filenames = sn_load_lines('specs_ms.txt');
end;
numSpecs = size(specs_raw,1);   sz=whos('specs_raw');

%
% Pepnovo filtering/scoring
%
if initial_stage <= STAGE_SCORING
	c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
	curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
    fprintf(1,'Starting stage: scoring (%s %s)\n',curDate,curTime);

    % Fake high charges (4+) as 3+ for Pepnovo scoring
    charges = zeros(numSpecs,1);   
    for s=1:numSpecs 
        if ~isempty(specs_raw{s,4}) 
            charges(s)=specs_raw{s,4};
            if specs_raw{s,4}>3 specs_raw{s,4}=3; end;
        end;
    end;
    sn_save_mgf(specs_raw(:,2), specs_raw(:,3), specs_raw(:,4), 'specs_ms.mgf');
    for s=1:numSpecs if charges(s)>3 specs_raw{s,4}=charges(s); end; end;
    
    % Pepnovo v2
%     if ~isempty(params.MIN_SPECTRUM_QUALITY) pepnovoMode = sprintf('-pmcsqs_and_prm %s',params.MIN_SPECTRUM_QUALITY); else pepnovoMode = '-prm_only'; end;
%     status = dos(sprintf('%s/PepNovo_bin -model_dir %s/Models -model LTQ_LOW_TRYP -file specs_ms.mgf -digest NON_SPECIFIC %s > specs_scored.prms',params.EXE_DIR,params.EXE_DIR,pepnovoMode)); if status~=0 fprintf(1,'ERROR executing %s/Pepnovo_bin!\n',params.EXE_DIR); res=-1; return; end;
% 	  status = dos(sprintf('%s/convert prms specs_scored.prms specs_scored.pklbin kept_indices.bin',params.EXE_DIR)); if status~=0 fprintf(1,'ERROR reading Pepnovo output!\n'); res=-1; return; end;

    % Pepnovo v3
    if ~isempty(params.MIN_SPECTRUM_QUALITY) pepnovoMode = sprintf('-min_filter_prob %s',params.MIN_SPECTRUM_QUALITY); else pepnovoMode = '-no_quality_filter'; end;
    if strcmp(lower(params.CORRECT_PM),'yes') pepnovoMode=strcat(pepnovoMode,' -correct_pm'); else pepnovoMode=strcat(pepnovoMode,' -use_spectrum_mz'); end;
    if strcmp(lower(params.GUESS_CHARGE),'no') pepnovoMode=strcat(pepnovoMode,' -use_spectrum_charge'); end;
    if ~isempty(params.PEPNOVO_PTMS) s=sprintf('%s -PTMs M+16:C+57:%s',pepnovoMode,params.PEPNOVO_PTMS); pepnovoMode=s; else pepnovoMode=strcat(pepnovoMode,' -PTMs M+16:C+57'); end;
    pepnovoCmd = sprintf('%s/PepNovo_bin -prm_norm -model_dir %s/Models_pepnovo -model %s -file specs_ms.mgf -fragment_tolerance %f -digest NON_SPECIFIC %s > specs_scored.prms\n',params.EXE_DIR,params.EXE_DIR,params.PEPNOVO_MODEL,peakTol,pepnovoMode);
    fprintf(1,pepnovoCmd);
    status = dos(pepnovoCmd); if status~=0 fprintf(1,'ERROR executing %s/Pepnovo_bin!\n',params.EXE_DIR); res=-1; return; end;
    status = dos(sprintf('%s/convert prmsv3 specs_scored.prms specs_scored_raw',params.EXE_DIR)); if status~=0 fprintf(1,'ERROR reading Pepnovo output!\n'); res=-1; return; end;
%     delete('specs_ms.mgf');
    
    % Adjustments to pepnovo's output
	specs = cell(numSpecs,5);
	idxKept = double(sn_load_binarray('specs_scored_raw.bin','uint32'))+1;
    specs(idxKept(:,1),:) = sn_load_pklbin('specs_scored_raw.pklbin',1);
    clear idxKept;
	for s=1:numSpecs 
        if ~isempty(specs{s,2}) specs{s,2}(:,2)=specs{s,2}(:,2)+1; specs{s,2}=specs{s,2}(find(specs{s,2}(:,2)>0),:); end; 
        if isempty(specs{s,3}) specs{s,3}=0; specs{s,4}=0; end;
        if specs_raw{s,4}>3 specs{s,4}=specs_raw{s,4}; end;   % Recover correct precursor charges
    end;
    sn_save_pklbin(specs, 'specs_scored.pklbin');
else
    specs = sn_load_pklbin('specs_scored.pklbin',1);
end;
% if sz.bytes>50000000 clear specs_raw; specs_raw=[];
%     fprintf(1,'Cleared specs_raw...\n'); 
% end;
% specs = sn_load_pklbin('specs_scored.pklbin',1);    specs_name = 'specs_scored';
% idxKept = double(sn_load_binarray('kept_indices.bin','uint32'))+1;
% specs_raw = specs_raw(idxKept,:);    filenames = filenames(idxKept,:);    clear idxKept;

%
% Alignment and spectral-stars consensus
%
cd('../aligns');
if initial_stage <= STAGE_FILTERPAIRS
	c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
	curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
    fprintf(1,'Starting stage: filterpairs (%s %s)\n',curDate,curTime);

    gridnodes = str2num(params.GRID_NUMNODES);
    stgParams = {sprintf('AA_DIFF_COUNT=%s',params.AA_DIFF_COUNT); sprintf('TOLERANCE_PEAK=%s',params.TOLERANCE_PEAK); sprintf('TOLERANCE_PM=%s',params.TOLERANCE_PM); sprintf('MAX_SHIFT=%s',params.MAX_MOD_MASS); sprintf('MIN_RATIO=%s',params.MIN_RATIO); 'MIN_SHIFT=0'; sprintf('MIN_NUM_MATCHED_PEAKS=%s',params.MIN_MATCHED_PEAKS)};
    if ~isempty(params.TAGS_FILTER) stgParams = [stgParams; {sprintf('INPUT_TAGS=%s',params.TAGS_FILTER);sprintf('TAGS_MATCH_FLANK=%s',params.TAGS_MATCH_FLANK);sprintf('TAGS_MATCH_COUNT=%s',params.TAGS_MATCH_COUNT)}]; end;
    if ~isempty(params.MIN_OVERLAP_AREA) stgParams = [stgParams; {sprintf('MIN_OVERLAP_AREA=%s',params.MIN_OVERLAP_AREA)}]; end;
    if strcmp(runMode,'sps') stgParams = [stgParams; {'PARTIAL_OVERLAPS=1'}]; end;
    if gridnodes<=0
        stgParams = [stgParams; {'INPUT_SPECS_PKLBIN=../spectra/specs_scored.pklbin'}];
        if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../%s',params.AMINO_ACID_MASSES)}]; end;
		fid = fopen('filterpairs.params','w');    if fid<=0 fprintf(1,'ERROR opening filterpairs.params!\n'); res=-1; return; end;
		fprintf(fid,'%d\nOUTPUT_ALIGNS=pairs_raw.bin\nMIN_PVALUE=%s\nOUTPUT_MEANS=means.bin\nOUTPUT_VARIANCE=vars.bin\nOUTPUT_RATIOS=ratios_raw.bin\n',5+size(stgParams,1),params.MIN_PVALUE);
		for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
		fclose(fid);
		status = dos(sprintf('%s/filterpairs filterpairs.params',params.EXE_DIR));   if status~=0 fprintf(1,'ERROR executing %s/filterpairs!\n',params.EXE_DIR); res=-1; return; end;

        fid = fopen('fn_aligns.txt','w');    if fid<=0 fprintf(1,'ERROR opening fn_aligns.txt!\n'); res=-1; return; end;
        fprintf(fid,'pairs_raw.bin\n');   fclose(fid);
        fid = fopen('fn_means.txt','w');     if fid<=0 fprintf(1,'ERROR opening fn_means.txt!\n'); res=-1; return; end;
        fprintf(fid,'means.bin\n');       fclose(fid);
        fid = fopen('fn_vars.txt','w');      if fid<=0 fprintf(1,'ERROR opening fn_vars.txt!\n'); res=-1; return; end;
        fprintf(fid,'vars.bin\n');        fclose(fid);
        fid = fopen('fn_ratios.txt','w');    if fid<=0 fprintf(1,'ERROR opening fn_ratios.txt!\n'); res=-1; return; end;
        fprintf(fid,'ratios_raw.bin\n');  fclose(fid);
        
        fid = fopen('mergefilter.params','w');    if fid<=0 fprintf(1,'ERROR opening mergefilter.params!\n'); res=-1; return; end;
        stgParams = {};
        if strcmp(params.FILTER_TRIGS,'yes') stgParams = [stgParams; {'FILTER_TRIGS=1'}]; end;
        if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../../%s',params.AMINO_ACID_MASSES)}]; end;
        if strcmp(runMode,'sps') stgParams = [stgParams; {'PARTIAL_OVERLAPS=1'}]; end;
        fprintf(fid,'%d\nINPUT_ALIGNS_MULTIPLE=fn_aligns.txt\nINPUT_MEANS_MULTIPLE=fn_means.txt\nINPUT_VARS_MULTIPLE=fn_vars.txt\nINPUT_RATIOS_MULTIPLE=fn_ratios.txt\n',10+size(stgParams,1));
        fprintf(fid,'MIN_PVALUE=%s\nOUTPUT_ALIGNS=pairs.bin\nOUTPUT_PVALUES=pvalues.bin\nOUTPUT_MEANS_STDDEVS=means-stddevs.bin\nOUTPUT_RATIOS=ratios.bin\nOUTPUT_INDICES=idxKeptPairs.bin\n',params.MIN_PVALUE);
        for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
        fclose(fid);
        
        status = dos(sprintf('%s/mergefilter mergefilter.params',params.EXE_DIR));   if status~=0 fprintf(1,'ERROR executing %s/mergefilter!\n',params.EXE_DIR); res=-1; return; end;
    else
        dirContents = dir('.'); dirContents = {dirContents(:).name}';
        stgParams = [stgParams; {'INPUT_SPECS_PKLBIN=specs_scored.pklbin'}];
        if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=%s',params.AMINO_ACID_MASSES)}]; end;
        if isempty(find(strcmp('grid',dirContents)))
            status = dos('mkdir grid');   cd('grid');
            fprintf(1,'The following files need to be copied to a single folder on your grid before executing run_jobs.sh\n - spectra/specs_scored.pklbin\n - all files in aligns/grid\n');
            if ~isempty(params.TAGS_FILTER) fprintf(1,' - %s\n',params.TAGS_FILTER); end;
            if ~isempty(params.AMINO_ACID_MASSES) fprintf(1,' - %s\n',params.AMINO_ACID_MASSES); end;
            sn_grid_prep(numSpecs, str2num(params.GRID_NUMNODES), str2num(params.GRID_NUMCPUS), params.GRID_EXE_DIR, stgParams);
            fprintf(1,'After the grid execution completes, resume the analysis with "%s(''%s'',''filterpairs'');" from your Matlab prompt or "%s %s filterpairs" from your command line (without " quotations marks)\n',runMode,paramsFN,runMode,paramsFN);
            cd('../..');
            res=0; return;
        else
            % Load the resulting files, compute/enforce p-values
            cd('grid');
		    fid = fopen('mergefilter.params','w');    if fid<=0 fprintf(1,'ERROR opening mergefilter.params!\n'); res=-1; return; end;
		    stgParams = {};
            if strcmp(params.FILTER_TRIGS,'yes') stgParams = [stgParams; {'FILTER_TRIGS=1'}]; end;
            if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../../%s',params.AMINO_ACID_MASSES)}]; end;
            if strcmp(runMode,'sps') stgParams = [stgParams; {'PARTIAL_OVERLAPS=1'}]; end;
            fprintf(fid,'%d\nINPUT_ALIGNS_MULTIPLE=fn_aligns.txt\nINPUT_MEANS_MULTIPLE=fn_means.txt\nINPUT_VARS_MULTIPLE=fn_vars.txt\nINPUT_RATIOS_MULTIPLE=fn_ratios.txt\n',10+size(stgParams,1));
		    fprintf(fid,'MIN_PVALUE=%s\nOUTPUT_ALIGNS=../pairs.bin\nOUTPUT_PVALUES=../pvalues.bin\nOUTPUT_MEANS_STDDEVS=means-stddevs.bin\nOUTPUT_RATIOS=../ratios.bin\nOUTPUT_INDICES=../idxKeptPairs.bin\n',params.MIN_PVALUE);
    		for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
    		fclose(fid);

            status = dos(sprintf('%s/mergefilter mergefilter.params',params.EXE_DIR));   if status~=0 fprintf(1,'ERROR executing %s/mergefilter!\n',params.EXE_DIR); res=-1; return; end;
            cd('..');
        end;
    end
end;
% pairs = sn_load_aligns('pairs.bin',pairsNumCols);

% specalign
if initial_stage <= STAGE_ALIGNMENT
	c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
	curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
    fprintf(1,'Starting stage: alignment (%s %s)\n',curDate,curTime);

    stgParams = {};
	if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../%s',params.AMINO_ACID_MASSES)}]; end;
    if strcmp(runMode,'sps') stgParams = [stgParams; {'PARTIAL_OVERLAPS=1'}]; end;
    fid = fopen('specalign.params','w');    if fid<=0 fprintf(1,'ERROR opening specalign.params!\n'); res=-1; return; end;
	fprintf(fid,'%d\nINPUT_SPECS_PKLBIN=../spectra/specs_scored.pklbin\nINPUT_ALIGNS=pairs.bin\nOUTPUT_SPECS=pairs.pklbin\nOUTPUT_MODPOS=pairs_modpos.bin\n',11+size(stgParams,1));
	fprintf(fid,'OUTPUT_STARS=../spectra/stars_only.pklbin\nOUTPUT_STARS_INDEX=../spectra/stars_indices.bin\nOUTPUT_STARS_ALL=../spectra/stars.pklbin\n');
	fprintf(fid,'TOLERANCE_PEAK=%s\nMAX_AA_JUMP=%s\nPENALTY_PTM=%s\nPENALTY_SAME_VERTEX=%s\n',params.TOLERANCE_PEAK,params.MAX_AA_JUMP,params.PENALTY_PTM,params.PENALTY_SAME_VERTEX);
    for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
	fclose(fid);
	
	status = dos(sprintf('%s/starden specalign.params',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/specalign!\n',params.EXE_DIR); res=-1; return; end;
end;
sz=whos('specs'); if sz.bytes>50000000 clear specs; 
    fprintf(1,'Cleared specs...\n'); 
end;
stars = sn_load_pklbin('../spectra/stars.pklbin');   stars_indices = sn_load_binarray('../spectra/stars_indices.bin','uint')+1;

% filterstarpairs
if initial_stage <= STAGE_FILTERSTARPAIRS
	c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
	curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
    fprintf(1,'Starting stage: filterstarpairs (%s %s)\n',curDate,curTime);

    % params file for filterstarpairs should include params.TOLERANCE_PM and params.MAX_MOD_MASS
    
    stgParams = {sprintf('MIN_MATCHED_PEAKS=%s',params.MIN_MATCHED_PEAKS)};
    if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../%s',params.AMINO_ACID_MASSES)}]; end;
    if strcmp(runMode,'sps') stgParams = [stgParams; {'PARTIAL_OVERLAPS=1'}]; end;
    fid = fopen('filterstarpairs.params','w');   if fid<=0 fprintf(1,'ERROR opening filterstarpairs.params!\n'); res=-1; return; end;
	fprintf(fid,'%d\nINPUT_SPECS_PKLBIN=../spectra/stars.pklbin\nINPUT_ALIGNS=pairs.bin\nOUTPUT_ALIGNS=pairs_stars.bin\n',12+size(stgParams,1));
	fprintf(fid,'TOLERANCE_PEAK=%s\nTOLERANCE_PM=%s\nMAX_MOD_MASS=%s\nMAX_AA_JUMP=%s\nPENALTY_PTM=%s\nPENALTY_SAME_VERTEX=%s\nMIN_RATIO=%s\n',params.TOLERANCE_PEAK,params.TOLERANCE_PM,params.MAX_MOD_MASS,params.MAX_AA_JUMP,params.PENALTY_PTM,params.PENALTY_SAME_VERTEX,params.MIN_RATIO);
	fprintf(fid,'OUTPUT_SPECS=none\nOUTPUT_MODPOS=none\n');
    for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
	fclose(fid);
	status = dos(sprintf('%s/filterstarpairs filterstarpairs.params',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/filterstarpairs!\n',params.EXE_DIR); res=-1; return; end;
    
    if strcmp(runMode,'sps')  % Remove networks containing only spectra from the same peptide (same parent masses and small shifts)
        pairs_stars = sn_load_aligns('pairs_stars.bin',pairsNumCols);
        idx = sps_filter_samepep([specs{:,3}]', pairs_stars, 2*pmTol+0.0001);
        pairs_stars = pairs_stars(idx,:);
        sn_save_aligns('pairs_stars.bin',pairs_stars(:,1:2),pairs_stars(:,3:6));
    end;   
end;
pairs_stars = sn_load_aligns('pairs_stars.bin',pairsNumCols);    cd('..');

%
% InsPecT
%
peptides = repmat({''},numSpecs,1);   proteins = repmat({''},numSpecs,1);
if ~isempty(params.INSPECT_PEPTIDES) | ~isempty(params.TEXT_PEPTIDES)
    if str2num(params.CLUSTER_MIN_SIZE)>0
        fprintf(1,'WARNING: InsPecT results from non-clustered spectra are not usable - InsPecT results should have been obtained by searching spectra/specs_ms\n');
    end;
    if ~isempty(params.INSPECT_PEPTIDES) 
        peptidesDBS = sn_load_inspect(peptides, params.INSPECT_PEPTIDES, filenames);
        proteinsDBS = peptidesDBS(:,6);    peptidesDBS = peptidesDBS(:,1);
    end;
    if ~isempty(params.TEXT_PEPTIDES) 
        peptidesDBS = sn_load_lines(params.TEXT_PEPTIDES); 
        if size(peptidesDBS,1)~=numSpecs fprintf(1,'ERROR: Number of lines in %s (%d) is different from number of spectra (%d)! Ignoring peptide annotations.\n',params.TEXT_PEPTIDES,size(peptidesDBS,1),numSpecs); peptidesDBS = peptides; end;
        proteinsDBS = proteins;
    end;
    save peptidesDBS peptidesDBS;
    save proteinsDBS proteinsDBS;
else peptidesDBS = peptides; proteinsDBS = proteins; end

%
% Shotgun Protein Sequencing
%
cd('assembly');
if strcmp(runMode,'sps')
    if initial_stage <= STAGE_ASSEMBLY
		c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
		curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
        fprintf(1,'Starting stage: assembly (%s %s)\n',curDate,curTime);
        
        stgParams = [{sprintf('MIN_MATCHED_PEAKS=%s',params.MIN_MATCHED_PEAKS)}; {sprintf('MIN_EDGES_TO_COMPONENT=%s',params.SPS_MIN_EDGES_TO_COMPONENT)}];
        stgParams = [stgParams; {sprintf('MAX_MOD_MASS=%s',params.MAX_MOD_MASS)}];
        if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../%s',params.AMINO_ACID_MASSES)}]; end;
        fid = fopen('assembly.params','w');   if fid<=0 fprintf(1,'ERROR opening assembly.params!\n'); res=-1; return; end;
		fprintf(fid,'%d\nINPUT_SPECS_PKLBIN=../spectra/stars.pklbin\nINPUT_ALIGNS=../aligns/pairs_stars.bin\n',11+size(stgParams,1));
		fprintf(fid,'TOLERANCE_PEAK=%s\nMAX_AA_JUMP=%s\nPENALTY_PTM=%s\nPENALTY_SAME_VERTEX=%s\nGRAPH_TYPE=2\n',params.TOLERANCE_PEAK,params.MAX_AA_JUMP,params.PENALTY_PTM,params.PENALTY_SAME_VERTEX);
		fprintf(fid,'OUTPUT_SPECS=sps_seqs.pklbin\nOUTPUT_MODPOS=sps_seqs_modpos.bin\nPATH_MIN_PEAKS=%s\nPATH_MIN_SPECS=%s\n', params.SPSPATH_MIN_NUM_PEAKS, params.SPSPATH_MIN_NUM_SPECS);
        for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
		fclose(fid);
		[status,output] = dos(sprintf('%s/masab assembly.params',params.EXE_DIR),'-echo');    if status~=0 fprintf(1,'ERROR executing %s/masab!\n',params.EXE_DIR); res=-1; return; end;
        sn_save_lines('../log_masab.txt',{output},0);
    end
    sps_seqs = sn_load_pklbin('sps_seqs.pklbin',1);
    sps_seqs_info  = sps_load_abinfo('component_info.bin');
    sps_seqs_comps = sps_load_bla('components.bla','int32');
end
cd('..');

if strcmp(runMode,'specnets')
	%
	% Tag generation + database search
	%
	cd('spectra');    minEI = str2num(params.MIN_PERC_EXPINT);    minTP = str2num(params.MIN_PERC_TP);
	if ~isempty(params.FASTA_DATABASE)
		if initial_stage <= STAGE_TAGSEARCH
			c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
			curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
            fprintf(1,'Starting stage: tagsearch (%s %s)\n',curDate,curTime);
		
            sn_save_lines('db_peptides.txt', peptidesDBS, 1);
            if str2num(params.TAG_LEN)>0
                stgParams = [{'SUBSET_IDX=stars_indices.bin'}; {'OUTPUT_MIDX=tagsearch_midx.pklbin'}];
                stgParams = [stgParams; {'OUTPUT_MATCHED_PROTS=tagsearch_mp.bin'}; {'OUTPUT_SPECS=tagsearch_outspecs.pklbin'}];
                stgParams = [stgParams; {'OUTPUT_EITP=tagsearch_eitp.bin'}; {'INPUT_PEPTIDES=db_peptides.txt'}];
                stgParams = [stgParams; {sprintf('MIN_PERC_EXPINT=%f',minEI)}; {sprintf('MIN_PERC_TP=%f',minTP)}];
                if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../%s',params.AMINO_ACID_MASSES)}]; end;
                fid = fopen('tagsearch.params','w');   if fid<=0 fprintf(1,'ERROR opening tagsearch.params!\n'); res=-1; return; end;
				fprintf(fid,'%d\nINPUT_SPECS_PKLBIN=stars.pklbin\nINPUT_FASTA=../%s\nOUTPUT_PEPTIDES=peptides_tagsearch.txt\n',10+size(stgParams,1),params.FASTA_DATABASE);
				fprintf(fid,'TOLERANCE_PEAK=%s\nTOLERANCE_PM=%s\nTAG_LEN=%s\nDOUBLE_AA_JUMPS=%s\nMATCH_TAG_FLANKING_MASSES=%s\nMAX_PARSIMONY=1\nMAX_NUM_TAGS=%s\n',params.TOLERANCE_PEAK,params.TOLERANCE_PM,params.TAG_LEN,params.DOUBLE_AA_JUMPS,params.MATCH_TAG_FLANKING_MASSES,params.MAX_NUM_TAGS);
                for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
				fclose(fid);
				status = dos(sprintf('%s/tagsearch tagsearch.params',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/tagsearch!\n',params.EXE_DIR); res=-1; return; end;
            else fid = fopen('peptides_tagsearch.txt','w'); fclose(fid); end;  % Just create empty file, use loaded peptide annotations (e.g. InsPecT)
%             delete('db_peptides.txt');
        end;
	end;
% 	eitp = double(sn_load_binarray('tagsearch_eitp.bin','uint16'))/10000; if size(stars,1)~=size(eitp,1) fprintf(1,'ERROR loading necessary data (eitp) - rerun with INITIAL_STAGE=tagsearch or earlier (in %s).\n',paramsFN); res=-1; return; end;
	tagsearch_peptides = sn_load_tagsearch(stars,'peptides_tagsearch.txt');
% 	tagsearch_peptides = {};
    stars(:,5) = peptidesDBS;   proteins = proteinsDBS;   % Initialize with database IDs
	if ~isempty(tagsearch_peptides)
		for pivot=1:length(tagsearch_peptides) 
            if tagsearch_peptides{pivot,2}==1 
                stars{pivot,5}=tagsearch_peptides{pivot,3}{1}; 
                proteins{pivot}=tagsearch_peptides{pivot,1}{1,5}; 
            end; 
        end;
		peptides = stars(:,5);
	end;
% 	stars_dbannot = sn_load_pklbin('../spectra/tagsearch_outspecs.pklbin');   stars_dbannot(:,5) = stars(:,5);
%     if size(stars_dbannot,1)~=size(peptides,1) fprintf(1,'ERROR loading necessary data (stars_dbannot) - rerun with INITIAL_STAGE=tagsearch or earlier (in %s).\n',paramsFN); res=-1; return; end;
	
	%
	% Spectral Networks: propagation of annotations
	%
	cd('../specnets');
    if initial_stage <= STAGE_SPECNETS
		c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
		curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
        fprintf(1,'Starting stage: specnets (%s %s)\n',curDate,curTime);
	
%         [stars_dbannot, eitp, proteins] = sn_enforceannots(stars, peptides, proteins, peptidesDBS, proteinsDBS, 'prm', peakTol, str2num(params.MAX_MOD_MASS), minEI, minTP);
% 		peptides_dbannot=stars_dbannot(:,5);   save peptides_dbannot peptides_dbannot;
% 		save proteins proteins;
        
        stgParams = [{'OUTPUT_SPECS_MATCHIDX=snets_midx.pklbin'}; {'OUTPUT_ALIGNS=snets_pairs.bin'}];
        if ~isempty(params.AMINO_ACID_MASSES) stgParams = [stgParams; {sprintf('AMINO_ACID_MASSES=../%s',params.AMINO_ACID_MASSES)}];  end;
        stgParams = [stgParams; {sprintf('TOLERANCE_PEAK=%s',params.TOLERANCE_PEAK)}];
		fid = fopen('pathproj.params','w');  if fid<=0 fprintf(1,'ERROR opening pathproj.params!\n'); res=-1; return; end;
		fprintf(fid,'%d\nINPUT_SPECS_PKLBIN=../spectra/tagsearch_outspecs.pklbin\nINPUT_ALIGNS=../aligns/pairs_stars.bin\nINPUT_ANNOTATED=../spectra/tagsearch_eitp.bin\nOUTPUT_ANNOTINFO=nets_annots.bin\n',8+size(stgParams,1));
		fprintf(fid,'MIN_PERC_EXPINT=%.0f\nMIN_PERC_TP=%.0f\nOUTPUT_SPECS_PROJ=snets_specs.pklbin\nMIN_MATCHED_PEAKS=%s\n',10000*minEI, 10000*minTP,params.MIN_MATCHED_PEAKS);
        for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
		fclose(fid);
		status = dos(sprintf('%s/pathproj pathproj.params',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/pathproj!\n',params.EXE_DIR); res=-1; return; end;
    end
	[stars_dbannot, pp_matches, pp_all] = sn_load_pathproj(stars, minEI, minTP, str2num(params.MAX_MOD_MASS), 'nets_annots.bin', 'stars');
    maxPropLvl = max(pp_matches(:,2));
    for lvl=1:maxPropLvl
        idxLvl = find(pp_matches(:,2)==lvl);
        for pivot=1:length(idxLvl)
            proteins{pp_matches(idxLvl(pivot),1)} = proteins{pp_matches(idxLvl(pivot),3)};
        end
    end
save tmp_stats stars_dbannot pp_matches pp_all;

	%
	% Spsplot statistics
	%
    c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
    curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
    fprintf(1,'Starting stage: spsplot statistics (%s %s)\n',curDate,curTime);
    % Save peptide IDs in InsPecT tab-delimited format
    lines = repmat({''},numSpecs+1,1);  curLine=2;
    for specIdx=1:numSpecs
        if ~isempty(stars_dbannot{specIdx,5})
            lines{curLine}=sprintf('\t%d\t*.%s.*\t%d\t%s',specIdx-1,stars_dbannot{specIdx,5},stars_dbannot{specIdx,4},proteins{specIdx});
            curLine=curLine+1;
        end
    end
    lines=lines(1:curLine-1,:);   sn_save_lines('snets_peptides.txt',lines,1);
    if ~isempty(params.SPSPLOT_OUTDIR)
        % Run spsplot to obtain ion statistics
        stgParams = [{sprintf('OUTDIR=%s/specnets',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_FASTA=%s/%s',projDir,params.FASTA_DATABASE)}];
        stgParams = [stgParams; {sprintf('FILE_MS=%s/spectra/specs_ms.pklbin',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_STARSINDEX=%s/spectra/stars_indices.bin',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_STARS=%s/specnets/snets_specs.pklbin',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_PRM=%s/spectra/specs_scored.pklbin',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_INSPECT=%s/specnets/snets_peptides.txt',projDir)}];
%         stgParams = [stgParams; {sprintf('FILE_PAIRS=%s/aligns/pairs.bin',projDir)}];
%         stgParams = [stgParams; {sprintf('FILE_MATCHPAIRS=%s/aligns/pairs.pklbin',projDir)}];
%         stgParams = [stgParams; {sprintf('FILE_STARSPAIRS=%s/aligns/pairs_stars.bin',projDir)}];
        fid = fopen('spsplot_stats.params','w');  if fid<=0 fprintf(1,'ERROR opening spsplot_stats.params!\n'); res=-1; return; end;
        fprintf(fid,'%d\n',size(stgParams,1));  for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
        fclose(fid);
        status = dos(sprintf('%s/spsplot --p spsplot_stats.params --stats',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/spsplot --stats!\n',params.EXE_DIR); res=-1; return; end;
    end

    %
	% SVM-based FDR
	%
% Filter out bad IDs 
%   - includes MS/MS statistics because it's not doable inside pathproj
%   - soon to be replaced with SVM-based filters
%     eitp=zeros(numSpecs,4);
%     stats = load('stats-ms.txt');      eitp(stats(:,1),1:2)=[stats(:,2) max(stats(:,4:5)')'];
%     stats = load('stats-stars.txt');   eitp(stats(:,1),3:4)=[stats(:,2) max(stats(:,4:5)')'];
%     idxClear = find(eitp(:,1)<minEI*10000 | eitp(:,2)<minTP*10000 | eitp(:,3)<minEI*10000 | eitp(:,4)<minTP*10000);
    peptides = sn_load_inspect(repmat({''},numSpecs,1), 'snets_peptides_fdr001.txt', []);
    idxClear = setdiff(find(strcmp(stars_dbannot(:,5),'')==0),find(strcmp(peptides(:,1),'')));
fprintf(1,' -- Removed %d annotations as non-significant\n',length(idxClear));
    tagsearch_midx = sn_load_pklbin('../spectra/tagsearch_midx.pklbin');   tagsearch_mp=sn_load_binarray('../spectra/tagsearch_mp.bin','int');
    for specIdx=1:length(idxClear) 
        stars_dbannot{idxClear(specIdx),5}='';   tagsearch_midx{idxClear(specIdx),2}=[];   tagsearch_mp(idxClear(specIdx),1)=-1;
    end;
    sn_save_lines('snets_peptidesF.txt',stars_dbannot(:,5),1);
    sn_save_pklbin(tagsearch_midx,'peptides_midx.pklbin');   sn_save_binarray('peptides_mp.bin',tagsearch_mp,'int');

    %
	% Homology glue (homglue): merge networks matched to overlapping protein locations
	%
    c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
    curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
    fprintf(1,'Starting stage: homglue (%s %s)\n',curDate,curTime);
    
	stgParams = [{'INPUT_SPECS_PKLBIN=snets_specs.pklbin'}];
    stgParams = [stgParams; {'INPUT_SPECNETS_ALIGNS=snets_pairs.bin'}; {'INPUT_SPECNETS_MATCHED=snets_midx.pklbin'}];
    stgParams = [stgParams; {'INPUT_SPECNETS_ALIGNS_ALL=../aligns/pairs_stars.bin'}];
%     stgParams = [stgParams; {'INPUT_MATCHED_PEAKS_IDX=../spectra/tagsearch_midx.pklbin'}; {'INPUT_MATCHED_PROTS=../spectra/tagsearch_mp.bin'}];
    stgParams = [stgParams; {'INPUT_MATCHED_PEAKS_IDX=peptides_midx.pklbin'}; {'INPUT_MATCHED_PROTS=peptides_mp.bin'}];
    stgParams = [stgParams; {'OUTPUT_CSV=homglue_matches.txt'}; {'OUTPUT_SPECS=homglue_main.pklbin'}];
    stgParams = [stgParams; {'OUTPUT_MATCHES_REF=homglue_ref'}; {'OUTPUT_MATCHES_CSPS=homglue_main'}];
    stgParams = [stgParams; {'SPEC_TYPE_MSMS=0'}; {'GRAPH_TYPE=2'}; {'MIN_CONTIG_SET=1'}];
    stgParams = [stgParams; {sprintf('INPUT_FASTA=%s/%s',projDir,params.FASTA_DATABASE)}];
    stgParams = [stgParams; {sprintf('TOLERANCE_PEAK=%s',params.TOLERANCE_PEAK)}];
    stgParams = [stgParams; {sprintf('TOLERANCE_PM=%s',params.TOLERANCE_PM)}];
    stgParams = [stgParams; {sprintf('RESOLUTION=%s',params.RESOLUTION)}];
    stgParams = [stgParams; {sprintf('PENALTY_PTM=%s',params.PENALTY_PTM)}];
    stgParams = [stgParams; {sprintf('PENALTY_SAME_VERTEX=%s',params.PENALTY_SAME_VERTEX)}];
    fid = fopen('homglue.params','w');  if fid<=0 fprintf(1,'ERROR opening homglue.params!\n'); res=-1; return; end;
    fprintf(fid,'%d\n',size(stgParams,1));  for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
    fclose(fid);
% return;
    [status,output] = dos(sprintf('%s/homglue homglue.params',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/homglue!\n',params.EXE_DIR); res=-1; return; end;
	sn_save_lines('../log_homglue.txt',{output},0);
return;

    %
	% SPSPlot: generate HTML report
	%
    if ~isempty(params.SPSPLOT_OUTDIR)
        c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
        curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
        fprintf(1,'Starting stage: spsplot (%s %s)\n',curDate,curTime);

        stgParams = [{sprintf('OUTDIR=%s',params.SPSPLOT_OUTDIR)}];
        stgParams = [stgParams; {sprintf('FILE_FASTA=%s/%s',projDir,params.FASTA_DATABASE)}];
        stgParams = [stgParams; {sprintf('FILE_MS=%s/spectra/specs_ms.pklbin',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_INDEX=%s/spectra/input_index.txt',projDir)}];
        if str2num(params.CLUSTER_MIN_SIZE)>0
            stgParams = [stgParams; {sprintf('FILE_CLUSTER=%s/spectra/out/clust/clusters_0_*.clust',projDir)}];
            stgParams = [stgParams; {sprintf('FILE_CLUSTERMS=%s/spectra/specs_ms_*.pklbin',projDir)}];
            stgParams = [stgParams; {sprintf('FILE_CLUSTERSCAN=%s/spectra/specs_ms_*.bin',projDir)}];
        end
        stgParams = [stgParams; {sprintf('FILE_STARS=%s/specnets/snets_specs.pklbin',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_COMP=%s/specnets/component_info.bin',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_SEQS=%s/specnets/homglue_main.pklbin',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_MP=%s/specnets/homglue_main_mp.bin',projDir)}];
        stgParams = [stgParams; {sprintf('FILE_MIDX=%s/specnets/homglue_main_midx.pklbin',projDir)}];
        fid = fopen('spsplot.params','w');  if fid<=0 fprintf(1,'ERROR opening spsplot.params!\n'); res=-1; return; end;
        fprintf(fid,'%d\n',size(stgParams,1));  for pIdx=1:size(stgParams,1) fprintf(fid,'%s\n',stgParams{pIdx}); end;
        fclose(fid);
        status = dos(sprintf('%s/spsplot --p spsplot.params',params.EXE_DIR));    if status~=0 fprintf(1,'ERROR executing %s/spsplot!\n',params.EXE_DIR); res=-1; return; end;
    end;
    
    cd('..');
end
return;

% Note: no more peptides_dbannot variable - exactly the same as peptides (tagsearch already reconciles with InsPecT)

%
% Report results
%
if initial_stage <= STAGE_REPORT
    c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
	curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
    fprintf(1,'Starting stage: report (%s %s)\n',curDate,curTime);

    if strcmp(runMode,'sps')
save tmp_sps sps_seqs sps_seqs_comps stars sps_seqs_info;
        sps_report(sps_seqs, sps_seqs_comps, str2num(params.TAG_LEN), stars, sps_seqs_info, peakTol, 'report/');
    else
        cd('specnets');
%         [stars_dbannot, pp_matches, pp_all] = sn_load_pathproj(stars_dbannot, minEI, minTP, str2num(params.MAX_MOD_MASS), 'nets_annots.bin', 'stars');
        [stars_dbannot, pp_matches, pp_all] = sn_load_pathproj(stars, minEI, minTP, str2num(params.MAX_MOD_MASS), 'nets_annots.bin', 'stars');
		idx = min(find(pp_matches(:,2)>0));
        if ~isempty(idx) for matchIdx=idx:size(pp_matches,1) proteins{pp_matches(matchIdx,1)}=proteins{pp_matches(matchIdx,3)}; end; end;
peptides_annot = stars_dbannot(:,5);   peptides = stars(:,5);
save tmp_peptides peptides tagsearch_peptides peptides_annot;

        % Floatation (repositioning) of modifications and replacement of sequence extension/deletion masses 
        %  with the corresponding amino acid sequence additions/deletions (respectively)
        if ~isempty(params.FASTA_DATABASE)
%             status = dos(sprintf('%s/msmod spectra/specs_ms.pklbin report_specnets.txt %s report_specnets_float.txt',params.EXE_DIR,params.FASTA_DATABASE));
%             if status~=0 fprintf(1,'ERROR executing %s/msmod!\n',params.EXE_DIR); res=-1; return; end;

%             peps_float = sn_load_lines('float.txt');
%             for i=1:length(peps_float)
%                 [s,p]=strtok(peps_float{i},char(9));   p=strtok(p,char(9));
%                 stars_dbannot{str2num(s)+1,5} = p;
%             end
        end
        
        % Recover previous unreliable peptide annotations
		for specIdx=1:size(stars_dbannot,1) if isempty(stars_dbannot{specIdx,5}) stars_dbannot{specIdx,5}=peptides{specIdx}; end; end;
        sn_save_lines('dbg_peptides_final.txt',stars_dbannot(:,5),1);
        
		% Compute % b/y-ions and explained intensity
		if ~isempty(specs_raw)
            sn_evalspecs([specs_raw(:,1:4) stars_dbannot(:,5)],'msms', peakTol,'',[],str2num(params.MAX_MOD_MASS),0,0,'specs_raw_stats',0,0);
            load specs_raw_stats_sn_evalspecs expInt;
            eiMSMS = sum(expInt(:,1:4)')';
		else eiMSMS = []; end;
		
        sn_evalspecs(stars_dbannot,'prm', peakTol,'',[],str2num(params.MAX_MOD_MASS),0,1,'stars_dbannot',1,1);
		load stars_dbannot_sn_evalspecs expInt bypercs blists;
		eiPRM = [sum(expInt(:,1:3)')' max(bypercs(:,1:4)')'];
		
        cysMass = AAmasses(min(find(AAletters=='C')));   if isempty(cysMass) cysMass=160.0306482; end;

%         % Select only pairs between annotated spectra
%         markNE = strcmp(stars_dbannot(:,5),'')~=1;
%         pairs_stars_annot = pairs_stars(find(min(markNE(pairs_stars(:,1:2)'))'==1),:);
%         % Select only pairs of spectra with overlapping IDs
%         peptides_nomod = stars_dbannot(:,5);
%         for i=1:size(peptides_nomod,1) s=peptides_nomod{i}; if ~isempty(s) s=s(find(s>='A' & s<'Z')); peptides_nomod{i}=s; end; end;
%         toKeep = zeros(size(pairs_stars_annot,1),1);
%         for i=1:size(toKeep) 
%             s1=peptides_nomod{pairs_stars_annot(i,1)};   len1 = length(s1);
%             s2=peptides_nomod{pairs_stars_annot(i,2)};   len2 = length(s2);
%             if len1>=len2
%                 if strncmp(s1,s2,len2) | strncmp(s1(1+len1-len2:len1),s2,len2) toKeep(i)=1; end;
%             else if strncmp(s1,s2,len1) | strncmp(s1,s2(1+len2-len1:len2),len1) toKeep(i)=1; end; end;
%         end
%         pairs_stars_annot = pairs_stars_annot(find(toKeep==1),:);
%         sn_save_aligns('pairs_stars_annot.bin',pairs_stars_annot(:,1:2),[pairs_stars_annot(:,3:5) zeros(size(pairs_stars_annot,1),1)]);

% specs     stars_dbannot       - annotated star spectra
% aligns    pairs_stars_annot   - pairs between annotated spectra (after propagation). Same as input to homglue?
%                                  used to define components, orient edges from smaller->larger parent masses and generate input files for graphviz
% pp_matches                    - pathproj matches: spectrum index (col.1) and propagation level (col.2)
% scores        eiPRM           - explained intensity / percent true b/y in PRM spectra
% threshExpInt  minEI           - minimum explained intensity threshold (to report an ID)
% threshTP      minTP           - minimum percent true b/y threshold (to report an ID)
% graphvizCmd   params.GRAPHVIZ_CMD - not used anymore
% labelpyDir    params.LABELPY_DIR - not used anymore
% cysMass                       - only needed for Label.py - not used anymore
% filenames                     - spectrum filenames to show on the report
% proteins                      - protein names for matched proteins
% specsMS   specs_raw           - only needed for Label.py - not used anymore
% eiMSMS    eiMSMS              - explained intensity / percent true b/y in MS/MS spectra 

% Pre 2010/04/26
%         cd('..');
%         nets = sn_report_specnets(stars_dbannot, pairs_stars_annot, pp_matches, eiPRM, minEI, minTP, params.GRAPHVIZ_CMD, params.LABELPY_DIR, cysMass, filenames, proteins, specs_raw, eiMSMS, 'report_specnets.txt');

        snets = sps_load_bla('components.bla','int32');
        specs_nets = zeros(numSpecs,1); for idxNet=1:length(snets) specs_nets(snets{netIdx})=netIdx; end;

        % Select only pairs between annotated spectra
        keep = zeros(size(pairs_stars,1),1);
        for pairIdx=1:size(pairs_stars,1)
            if specs_snets(pairs_stars(pairIdx,1))==specs_snets(pairs_stars(pairIdx,2)) keep(pairIdx)=1; end;
        end
        pairs_stars_annot = pairs_stars(find(keep==1),:);
        sn_save_aligns('pairs_stars_annot.bin',pairs_stars_annot(:,1:2),[pairs_stars_annot(:,3:5) zeros(size(pairs_stars_annot,1),1)]);

        cd('..');
        sn_report_specnets(stars_dbannot, pairs_stars_annot, snets, pp_matches, eiPRM, minEI, minTP, params.GRAPHVIZ_CMD, params.LABELPY_DIR, cysMass, filenames, proteins, specs_raw, eiMSMS, 'report_specnets.txt');

%         c=clock;   curDate=sprintf('%4.0f/%2.0f/%2.0f',c(1),c(2),c(3));    curDate(find(curDate==' '))='0';
% 		curTime=sprintf('%2.0f:%2.0f:%2.0f',c(4),c(5),c(6));    curTime(find(curTime==' '))='0';
%         fprintf(1,'Starting stage: report/ABruijn (%s %s)\n',curDate,curTime);
% 
%         save tmp1 stars_dbannot nets blists;
% %         load tmp1;
% %         
% %         nets_idx = sort(nets{74});   % Non-floated sag_p1
% %         nets_idx = sort(nets{480});  % Floated sag_p1
% %         nets_idx = sort(nets{237});  % Non-floated sag_p2
% %         nets_idx = sort(nets{275});   % Non-floated sag_p2a
%         nets_idx = sort(nets{64});   % Non-floated sag_exp1c
%         nets_sz = length(nets_idx);   minNumPeaks = str2num(params.SPSPATH_MIN_NUM_PEAKS);
%         nets_contigs = repmat({[]},numSpecs,1);   nets_seqs = repmat({[],[],0,0,[]},numSpecs,1);
%         for cIdx=1:nets_sz
%             sIdx = nets_idx(cIdx);
%             nets_contigs{sIdx} = sIdx;
%             if size(blists{sIdx,2},1)>=minNumPeaks 
% %                 nets_seqs{sIdx,2} = stars_dbannot{sIdx,2}([blists{sIdx,2}{:}],:);
%                 nets_seqs{sIdx,2} = [cumsum(sn_getmasses(stars_dbannot{sIdx,5},'',[],str2num(params.MAX_MOD_MASS)))'];
%                 nets_seqs{sIdx,2} = [nets_seqs{sIdx,2} repmat(1,size(nets_seqs{sIdx,2},1),1)];
%             end;
%         end
% 
%         [nets_abinfo, nets_seqs] = sn_asAbruijn(stars_dbannot, nets_seqs, nets_contigs, blists, minNumPeaks);
%         sn_save_abinfo(nets_abinfo,'specnets/nets_abinfo.bin');
%         save_pklbin(nets_seqs,'specnets/nets_seqs.pklbin');
%         save tmp2 nets_contigs nets_seqs nets_abinfo;
    end
end;

res=0;
