function sn_save_aligns(filename,pairs,aligns)
% function sn_save_aligns(filename,pairs,aligns)

fid = fopen(filename,'w');
if fid<=0 fprintf(1,'Error opening %s!\n',filename); return; end

szAligns = size(pairs,1);   count = fwrite(fid,szAligns,'int32');   if count~=1 fprintf(1,'ERROR writing count to %s!\n',filename); fclose(fid); return; end;
numWritten = 0;
while numWritten<szAligns  % Output 100,000 spectral pairs at a time
    toWrite = min(100000,szAligns-numWritten);
    data = [single(pairs(numWritten+1:numWritten+toWrite,:)) single(aligns(numWritten+1:numWritten+toWrite,:))]; % Write zero-based spectrum indices
    count = fwrite(fid,data','float');   
    if count~=toWrite*size(data,2) fprintf(1,'ERROR writing spectral pairs to %s!\n',filename); fclose(fid); return; end;
    numWritten = numWritten+toWrite;
end
fclose(fid);
