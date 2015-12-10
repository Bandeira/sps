function save_resultsPA(filename,separator,pairs,aligns)
% function save_resultsPA(filename,separator,pairs,aligns)
%
%  Saves a set of PA (partial overlap) aligns for further cpplib processing (e.g. masab)
%

fid = fopen(filename,'w');
if fid<=0 fprintf(1,'Error opening %s!\n',filename); return; end

szAligns = size(pairs,1);  if size(aligns,2)>3 idxLast=4; else idxLast=3; end;
fprintf(fid,'%.0f\n',szAligns);
for i=1:szAligns
    fprintf(fid,'%.0f%s%.0f%s%.1f%s%.1f%s%.1f%s%.1f\n',double(pairs(i,1)),separator,double(pairs(i,2)),separator,double(aligns(i,1)),separator,double(aligns(i,2)),separator,double(aligns(i,3)),separator,double(aligns(i,idxLast)));
end;

fclose(fid);
