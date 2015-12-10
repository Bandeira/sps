function save_clusters_specset(specSet, filename, filterEndPeaks, maximizeEndpointScore)
% function save_clusters_specset(specSet, filename, filterEndPeaks, maximizeEndpointScore)
%
%  Saves a set of spectra as a clusters file (suitable for matchma)
%
%  filterEndPeaks - if ==1 then all peaks <=50Da or >=PM-50Da are removed from the spectrum
%

fid=fopen(filename,'w');  if fid<=0 fprintf(1,'ERROR opening %s!\n',filename); return; end;
if size(specSet,2)==5 idxSpec=2; idxPM=3; else idxSpec=3; idxPM=5; end;

numSpecs = size(specSet,1);
fprintf(fid,'%.0d\n',numSpecs);
for i=1:size(specSet,1)
    spec = specSet{i,idxSpec}; if filterEndPeaks idx=find(spec(:,1)>50 & spec(:,1)<specSet{i,idxPM}-50); spec=spec(idx,:); end;
    numPeaks = size(spec,1);   if maximizeEndpointScore maxScore=max(spec(:,2)); else maxScore=0; end;

    fprintf(fid,'0\n');  % No shifts
    fprintf(fid,'%d %.1f\n',numPeaks,specSet{i,idxPM});  % Num. peaks ; Parent mass
    for p=1:numPeaks fprintf(fid,'%.1f %.1f\n',spec(p,1:2)); end;
%     fprintf(fid,'0\n');  % No additional endpoints
    fprintf(fid,'2\n0 %.1f\n%.1f %.1f\n',maxScore,specSet{i,idxPM}-19,maxScore);  % Zero/PM-19 b-ion endpoints
end

fclose(fid);
