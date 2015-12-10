function sn_save_pklbin(specs, filename)
% function sn_save_pklbin(specs, filename)

fid=fopen(filename,'w');  if fid<=0 fprintf(1,'Error opening file %s!\n',filename); return; end;
if size(specs,2)==5 idxSpec=2; idxPM=3; idxZ=4; else idxSpec=3; idxPM=5; idxZ=-1; end;

numSpecs = size(specs,1);   sizes=zeros(numSpecs,1);
for s=1:numSpecs sizes(s)=size(specs{s,idxSpec},1); end;

fwrite(fid,numSpecs,'int32');   fwrite(fid,sizes,'int16');
for s=1:numSpecs
    if isempty(specs{s,idxSpec}) data = [specs{s,idxPM} 0]; else data = [specs{s,idxPM} 0; specs{s,idxSpec}(:,1:2)]; end;
    if idxZ>0 & ~isempty(specs{s,idxZ}) data(1,2)=specs{s,idxZ}; end;
    fwrite(fid,data','float32');
end

fclose(fid);
