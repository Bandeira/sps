function gpf_batch_clst(baseName, specsType, peakTol, pmTol, resolution, minRatio,paramsFN)
% function gpf_batch_clst(baseName, specsType, peakTol, pmTol, resolution, minRatio,paramsFN)
%
%  Generates a .params file for batch_clst.
%

fid=fopen(paramsFN,'w');   if fid<=0 fprintf(1,'ERROR opening %s\n',paramsFN); return; end;

fprintf(fid,'8\n');
if strcmp(specsType,'pkl') fprintf(fid,'INPUT_SPECS=%s.pkl\n',baseName); end;
if strcmp(specsType,'pklbin') fprintf(fid,'INPUT_SPECS_PKLBIN=%s.pklbin\n',baseName); end;

fprintf(fid,'OUTPUT_PAIRS=%s_pairs.bin\n',baseName);
fprintf(fid,'OUTPUT_RATIOS=%s_ratios.bin\n',baseName);
fprintf(fid,'OUTPUT_CLUSTERS=%s_clusters.bla\n',baseName);
fprintf(fid,'TOLERANCE_PEAK=%f\nTOLERANCE_PM=%f\nRESOLUTION=%f\nMIN_RATIO=%f\n',peakTol, pmTol, resolution, minRatio);

fclose(fid);
