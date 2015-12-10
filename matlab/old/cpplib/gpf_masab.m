function gpf_masab(paramsFN, specsFN, alignsFN, alignsPAFN, labelsFN, outputBaseFN, minRatio, masabParams)
% function gpf_masab(paramsFN, specsFN, alignsFN, alignsPAFN, labelsFN, outputBaseFN, minRatio, masabParams)
%
%  Generates masab parameter file paramsFN.
%
%  Mandatory params: specsFN, alignsFN, outputBaseFN, minRatio
%  Optional params:  alignsPAFN, labelsFN, masabParams
%
%  masabParams - column cell vector with values for the execution parameters:
%                 MAX_AA_JUMP=1, TOLERANCE_PEAK=0.5, PENALTY_PTM=-200, PENALTY_SAME_VERTEX=-1000000
%                 GRAPH_TYPE=2 (these default values are used when masabParams is empty)

fid = fopen(paramsFN,'w'); if fid<=0 fprintf(1,'Error opening parameters file: %s\n',paramsFN); return; end;

if isempty(masabParams)
    masabParams={'MAX_AA_JUMP=1';'TOLERANCE_PEAK=0.5';'PENALTY_PTM=-200';'PENALTY_SAME_VERTEX=-1000000';'GRAPH_TYPE=2'};
end

% Header
% fprintf(fid,'%.0d\n',4 + size(masabParams,1) + ~isempty(alignsPAFN) + ~isempty(labelsFN) + (minRatio>0)*(2+~isempty(alignsPAFN)));  % 4 is for 2 input + 2 output files
fprintf(fid,'%.0d\n',4 + size(masabParams,1) + ~isempty(alignsPAFN) + ~isempty(labelsFN) + 2+~isempty(alignsPAFN) );  % 4 is for 2 input + 2 output files

% Input files
idx=max(find(specsFN=='.'));   extension = specsFN(idx+1:length(specsFN));
if ~isempty(extension) & strcmp(extension,'pklbin') fprintf(fid,'INPUT_SPECS_PKLBIN=%s\n',specsFN);
else fprintf(fid,'INPUT_SPECS=%s\n',specsFN); end;
fprintf(fid,'INPUT_ALIGNS=%s\n',alignsFN);
if ~isempty(alignsPAFN) fprintf(fid,'INPUT_ALIGNSPA=%s\n',alignsPAFN); end;
if ~isempty(labelsFN) fprintf(fid,'INPUT_LABELS=%s\n',labelsFN); end;

% Output files
fprintf(fid,'OUTPUT_SPECS=%s_masab_specs.pklbin\nOUTPUT_MODPOS=%s_masab_modPos.txt\n',outputBaseFN, outputBaseFN);
fprintf(fid,'MIN_RATIO=%f\n',minRatio); 
fprintf(fid,'OUTPUT_RATIOS=%s_masab_ratios.bin\n',outputBaseFN); 
if ~isempty(alignsPAFN) fprintf(fid,'OUTPUT_RATIOSPA=%s_masab_ratiospa.bin\n',outputBaseFN); end;

% Params
for i=1:size(masabParams,1) fprintf(fid,'%s\n',masabParams{i}); end;

fclose(fid);
