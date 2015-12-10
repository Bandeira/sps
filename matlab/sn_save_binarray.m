function res = sn_save_binarray(filename, data, type)
% function res = sn_save_binarray(filename, data, type)

fid=fopen(filename,'w');  if fid<=0 fprintf(1,'Error opening file %s!\n',filename); res=-1; return; end;

fwrite(fid,size(data),'int32');
fwrite(fid,data',type);

fclose(fid);
res=0;
