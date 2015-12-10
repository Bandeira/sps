function gpf_batch_asppa(paramsFN, specsFN, outputBaseFN, otherParams)
% function gpf_batch_asppa(paramsFN, specsFN, outputBaseFN, otherParams)
%
%  Generates batch_asp parameter file paramsFN.
%
%  Mandatory params: specsFN, alignsFN, outputBaseFN
%  Optional params:  otherParams
%
%  otherParams - column cell vector with values for the execution parameters:
%                 {'AA_DIFF_COUNT=-1';'TOLERANCE_PEAK=0.5';'TOLERANCE_PM=1.0';'MIN_SHIFT=0';'MAX_SHIFT=100';'MIN_RATIO=0.4'} - for mods (default)
%                 {'AA_DIFF_COUNT=2';'TOLERANCE_PEAK=0.5';'TOLERANCE_PM=1.0';'MIN_SHIFT=0';'MIN_RATIO=0.4'} - for asp aligns
%                 {'IDX_START=0';'IDX_END=9999';'TOLERANCE_PEAK=0.5';'MIN_OVERLAP_AREA=0.5';'MIN_RATIO=0.4'} - for pa aligns
%

fid = fopen(paramsFN,'w'); if fid<=0 fprintf(1,'Error opening parameters file: %s\n',paramsFN); return; end;

if isempty(otherParams)
    otherParams={'AA_DIFF_COUNT=-1';'TOLERANCE_PEAK=0.5';'TOLERANCE_PM=1.0';'MIN_SHIFT=0';'MAX_SHIFT=100';'MIN_RATIO=0.4'};  % Default is to generate params file to search for mods of mass <=100 Da
end

% Header
fprintf(fid,'%.0d\n',5 + size(otherParams,1));  % 5 is for 1 input + 4 output files

% Input/Output files
fprintf(fid,'INPUT_SPECS_PKLBIN=%s\nOUTPUT_ALIGNS=%s_aligns.txt\nOUTPUT_RATIOS=%s_ratios.txt\nOUTPUT_MEANS=%s_means.txt\nOUTPUT_VARIANCE=%s_vars.txt\n',specsFN,outputBaseFN,outputBaseFN,outputBaseFN,outputBaseFN);

% Params
for i=1:size(otherParams,1) fprintf(fid,'%s\n',otherParams{i}); end;

fclose(fid);
