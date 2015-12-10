function [res, params] = sn_read_params_file(filename, setDefaults)
% function [res, params] = sn_read_params_file(filename, setDefaults)
%
%  Reads a parameters file of the format
%
%    PARAM_NAME=<param_value>
%
%  And outputs a struct params with fields params.PARAM_NAME=<param_value> .
%  Lines starting with '#' are ignored (in the params file).
%
%  setDefaults - if omitted or different from 0, default fields and values are added to params
%  res = 1 if file loaded ok, 0 otherwise
%

params = struct('TOLERANCE_PEAK','0.5','TOLERANCE_PM','1','RESOLUTION','0.1');
fid = fopen(filename,'r'); if fid<=0 fprintf(1,'ERROR: Could not open %s!\n',filename); res=0; return; end; res=1;

if nargin<2 setDefaults=1; end;
if setDefaults~=0
	params.INPUT_SPECS='';
	params.FASTA_DATABASE='';
	params.AMINO_ACID_MASSES='';
	params.GRAPHVIZ_CMD='';
	params.INITIAL_STAGE='';
	params.LABELPY_DIR = '';
    
	% Pepnovo/MS-Cluster/Inspect parameters
	params.CLUSTER_MIN_SIZE='0';
	params.CLUSTER_PMTOL_PPM='';
	params.PEPNOVO_MODEL='CID_IT_TRYP';
	params.PEPNOVO_PTMS='';
	params.MIN_SPECTRUM_QUALITY='0.15';
	params.CORRECT_PM='yes';
	params.GUESS_CHARGE='yes';
	
	% Third party peptide annotations
	params.TEXT_PEPTIDES='';
	params.INSPECT_PEPTIDES='';
	params.INSPECT_ENFORCE='no';
	
	% Comparative Shotgun Protein Sequencing (CSPS) parameters
	params.SPS_PROJECTS='';
	params.CLUSTALW_EXE_DIR='';
	params.CLUSTALW_MINSCORE=250;
	params.REPORT_MIN_TAG=6;
	params.FORCE_REFERENCE='-1';
	
	% Alignment parameters
	params.AA_DIFF_COUNT='2';
	params.MIN_SHIFT='0';
	params.MIN_MOD_MASS='-100';
	params.MAX_MOD_MASS='100';
    params.MAX_NUM_MODS='1';
	params.MIN_RATIO='0.4';
	params.MIN_PVALUE='0.05';
	params.MIN_MATCHED_PEAKS='4';
	params.MAX_AA_JUMP='1';
	params.MIN_OVERLAP_AREA='';
	params.PENALTY_PTM='-200';
	params.PENALTY_SAME_VERTEX='-1000000';
	params.FILTER_TRIGS='no';
	params.TAGS_FILTER='';   % Used to filter possible spectral pairs
	params.TAGS_MATCH_FLANK='1';
	params.TAGS_MATCH_COUNT='2';
	
	% Sequencing parameters
	params.DENOVO='';
	params.DENOVO_EVAL='';
    params.SPSPATH_MIN_NUM_PEAKS='5';
    params.SPSPATH_MIN_NUM_SPECS='2';
    params.SPS_MIN_EDGES_TO_COMPONENT='1';
	
	% tagsearch/matchma parameters
	params.TAG_LEN='6';
	params.DOUBLE_AA_JUMPS='1';
	params.MATCH_TAG_FLANKING_MASSES='2';
	params.MAX_NUM_TAGS='0';
    params.MIN_MATCHED_PEAKS_DB='6';
	
	% Networks parameters (pathproj)
	params.MIN_PERC_EXPINT='0.01';
	params.MIN_PERC_TP='0.01';
	
	% Grid parameters
	params.GRID_NUMNODES='0';
	params.GRID_NUMCPUS='2';
	params.GRID_EXE_DIR='';
    
    % Reporting parameters
    params.SPSPLOT_OUTDIR='';
end;

while 1
    line = fgetl(fid);
    if ~ischar(line) break; end;
    if isempty(line) | line(1)=='#' | ~isempty(str2num(line(1))) continue; end;
    
    [token,value] = strtok(line,'=');
    params = setfield(params,token,value(2:length(value)));
end

fclose(fid);

if aux_check_type_range('TOLERANCE_PEAK', params.TOLERANCE_PEAK, 0, 1)==-1 res=-1; return; end;
if aux_check_type_range('TOLERANCE_PM', params.MIN_SPECTRUM_QUALITY, 0, 3)==-1 res=-1; return; end;
if aux_check_type_range('RESOLUTION', params.MIN_SPECTRUM_QUALITY, 0.01, 1)==-1 res=-1; return; end;
if aux_check_type_range('MIN_SPECTRUM_QUALITY', params.MIN_SPECTRUM_QUALITY, 0, 1)==-1 res=-1; return; end;
if aux_check_type_range('AA_DIFF_COUNT', params.AA_DIFF_COUNT, -inf, inf)==-1 res=-1; return; end;
if aux_check_type_range('MIN_SHIFT', params.MIN_SHIFT, 0, inf)==-1 res=-1; return; end;
if aux_check_type_range('MAX_MOD_MASS', params.MAX_MOD_MASS, 0, inf)==-1 res=-1; return; end;
if aux_check_type_range('MAX_NUM_MODS', params.MAX_NUM_MODS, 0, inf)==-1 res=-1; return; end;
if aux_check_type_range('MIN_RATIO', params.MIN_RATIO, 0, 1)==-1 res=-1; return; end;
if aux_check_type_range('MIN_PVALUE', params.MIN_PVALUE, 0, 1)==-1 res=-1; return; end;
if aux_check_type_range('MIN_MATCHED_PEAKS', params.MIN_MATCHED_PEAKS, 1, inf)==-1 res=-1; return; end;
if aux_check_type_range('MAX_AA_JUMP', params.MAX_AA_JUMP, 0, 3)==-1 res=-1; return; end;
if aux_check_type_range('PENALTY_PTM', params.PENALTY_PTM, -inf, inf)==-1 res=-1; return; end;
if aux_check_type_range('PENALTY_SAME_VERTEX', params.PENALTY_SAME_VERTEX, -inf, inf)==-1 res=-1; return; end;
if aux_check_type_range('TAG_LEN', params.TAG_LEN, 1, inf)==-1 res=-1; return; end;
if aux_check_type_range('DOUBLE_AA_JUMPS', params.DOUBLE_AA_JUMPS, 0, 1)==-1 res=-1; return; end;
if aux_check_type_range('MATCH_TAG_FLANKING_MASSES', params.MATCH_TAG_FLANKING_MASSES, 0, 2)==-1 res=-1; return; end;
if aux_check_type_range('MIN_PERC_EXPINT', params.MIN_PERC_EXPINT, 0, 1)==-1 res=-1; return; end;
if aux_check_type_range('MIN_PERC_TP', params.MIN_PERC_TP, 0, 1)==-1 res=-1; return; end;


function status = aux_check_type_range(param_name, param_str, min_val, max_val)
% Check type and range of numeric input parameters

status=-1;
if isempty(str2num(param_str)) fprintf(1,'ERROR: %s must be a number!\n',param_name);  return; end;
value = str2num(param_str);
if value<min_val | value>max_val fprintf(1,'ERROR: %s is outside the allowed range [%f,%f]!\n',param_name,min_val,max_val);  return; end;
