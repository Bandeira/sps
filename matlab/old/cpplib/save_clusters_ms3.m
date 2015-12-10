function save_clusters_ms3(ms3sets, filename)
% function save_clusters_ms3(ms3sets, filename)
%
%  Saves a set of spectra as a clusters file (suitable for ms3den)
%

fid=fopen(filename,'w');  if fid<=0 fprintf(1,'ERROR opening %s!\n',filename); return; end;

numSets = size(ms3sets,1);
fprintf(fid,'%.0d\n',numSets);
for i=1:numSets
    fprintf(fid,'%d\n',size(ms3sets{i},2));  % Num. spectra
    fprintf(fid,'%d 0 0\n',ms3sets{i}-1);
    
    fprintf(fid,'0 0\n');  % Num. peaks ; Parent mass
    fprintf(fid,'0\n');  % No additional endpoints
end

fclose(fid);
