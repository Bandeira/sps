function gpf_pathproj(paramsFN, specsFN, alignsFN, annotsFN, minPercExpInt, minPercTP, outputFN)
% function gpf_pathproj(paramsFN, specsFN, alignsFN, annotsFN, outputFN)
%
%  Generates parameters file paramsFN for pathproj
%
%  INPUT_SPECS=specsFN
%  INPUT_ALIGNS=alignsFN
%  INPUT_ANNOTATED=annotsFN
%  MIN_PERC_EXPINT=minPercExpInt
%  MIN_PERC_TP=minPercTP
%  OUTPUT_ANNOTINFO=outputFN
%

fid = fopen(paramsFN,'w'); if fid<=0 fprintf(1,'Error opening parameters file: %s\n',paramsFN); return; end;
fprintf(fid,'6\nINPUT_SPECS=%s\nINPUT_ALIGNS=%s\nINPUT_ANNOTATED=%s\nOUTPUT_ANNOTINFO=%s\n',specsFN, alignsFN, annotsFN, outputFN);
fprintf(fid,'MIN_PERC_EXPINT=%.0d\nMIN_PERC_TP=%.0d\n',10000*minPercExpInt, 10000*minPercTP);
fclose(fid);
