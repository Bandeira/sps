function [res, params] = read_params_file(filename)
% function params = read_params_file(filename)
%
%  Reads a parameters file of the format
%
%    PARAM_NAME=<param_value>
%
%  And outputs a struct params with fields params.PARAM_NAME=<param_value> .
%  Lines starting with '#' are ignored (in the params file).
%
%  res = 1 if file loaded ok, 0 otherwise
%

params = struct('TOLERANCE_PEAK','0.5','TOLERANCE_PM','1','RESOLUTION','0.1');
fid = fopen(filename,'r'); if fid<=0 fprintf(1,'ERROR: Could not open %s!\n',filename); res=0; return; end; res=1;

params.INPUT_SPECS='';
params.GRAPHVIZ_CMD='';
params.INITIAL_STAGE='';

% Alignment parameters
params.MIN_SPECTRUM_QUALITY='0.15';
params.AA_DIFF_COUNT='2';
params.MIN_SHIFT='0';
params.MAX_MOD_MASS='100';
params.MIN_RATIO='0.4';
params.MIN_PVALUE='0.05';
params.MIN_MATCHED_PEAKS='4';
params.MAX_AA_JUMP='1';
params.PENALTY_PTM='-200';
params.PENALTY_SAME_VERTEX='-1000000';

% tagsearch parameters
params.TAG_LEN='6';
params.DOUBLE_AA_JUMPS='1';
params.MATCH_TAG_FLANKING_MASSES='2';

% Networks parameters (pathproj)
params.MIN_PERC_EXPINT='0';
params.MIN_PERC_TP='0';

while 1
    line = fgetl(fid);
    if ~ischar(line) break; end;
    if isempty(line) | line(1)=='#' continue; end;
    
    [token,value] = strtok(line,'=');
    params = setfield(params,token,value(2:length(value)));
end

fclose(fid);

if aux_check_type_range('TOLERANCE_PEAK', params.TOLERANCE_PEAK, 0, 1)==-1 status=-1; return; end;
if aux_check_type_range('TOLERANCE_PM', params.MIN_SPECTRUM_QUALITY, 0, 3)==-1 status=-1; return; end;
if aux_check_type_range('RESOLUTION', params.MIN_SPECTRUM_QUALITY, 0.01, 1)==-1 status=-1; return; end;
if aux_check_type_range('MIN_SPECTRUM_QUALITY', params.MIN_SPECTRUM_QUALITY, 0, 1)==-1 status=-1; return; end;
if aux_check_type_range('AA_DIFF_COUNT', params.AA_DIFF_COUNT, -inf, inf)==-1 status=-1; return; end;
if aux_check_type_range('MIN_SHIFT', params.MIN_SHIFT, 0, inf)==-1 status=-1; return; end;
if aux_check_type_range('MAX_MOD_MASS', params.MAX_MOD_MASS, 0, inf)==-1 status=-1; return; end;
if aux_check_type_range('MIN_RATIO', params.MIN_RATIO, 0, 1)==-1 status=-1; return; end;
if aux_check_type_range('MIN_PVALUE', params.MIN_PVALUE, 0, 1)==-1 status=-1; return; end;
if aux_check_type_range('MIN_MATCHED_PEAKS', params.MIN_MATCHED_PEAKS, 1, inf)==-1 status=-1; return; end;
if aux_check_type_range('MAX_AA_JUMP', params.MAX_AA_JUMP, 0, 3)==-1 status=-1; return; end;
if aux_check_type_range('PENALTY_PTM', params.PENALTY_PTM, -inf, inf)==-1 status=-1; return; end;
if aux_check_type_range('PENALTY_SAME_VERTEX', params.PENALTY_SAME_VERTEX, -inf, inf)==-1 status=-1; return; end;
if aux_check_type_range('TAG_LEN', params.TAG_LEN, 1, inf)==-1 status=-1; return; end;
if aux_check_type_range('DOUBLE_AA_JUMPS', params.DOUBLE_AA_JUMPS, 0, 1)==-1 status=-1; return; end;
if aux_check_type_range('MATCH_TAG_FLANKING_MASSES', params.MATCH_TAG_FLANKING_MASSES, 0, 2)==-1 status=-1; return; end;
if aux_check_type_range('MIN_PERC_EXPINT', params.MIN_PERC_EXPINT, 0, 1)==-1 status=-1; return; end;
if aux_check_type_range('MIN_PERC_TP', params.MIN_PERC_TP, 0, 1)==-1 status=-1; return; end;


function status = aux_check_type_range(param_name, param_str, min_val, max_val)
% Check type and range of numeric input parameters

status=-1;
if isempty(str2num(param_str)) fprintf(1,'ERROR: %s must be a number!\n',param_name);  return; end;
value = str2num(param_str);
if value<min_val | value>max_val fprintf(1,'ERROR: %s is out of range!\n',param_name);  return; end;
