function sn_save_mgf(specs, filename, filenames)
% function sn_save_mgf(specs, parentMasses, parentCharges, filename, scans, filenames, peptides)

if nargin<3 filenames={}; end;

fid = fopen(filename,'w'); if fid<0 fprintf(1,'Error opening %s!\n',filename); end;
numSpecs = size(specs,1);
for i=1:numSpecs
    if isempty(specs{i,1}) scan = '';
    else scan = sprintf('SCANS=%.0d',specs{i,1}); end

    if isempty(specs{i,5}) seq = '';
    else seq = sprintf('SEQ=%s',specs{i,5}); end

    if isempty(filenames) | i>size(filenames,1) title = '';
    else title = sprintf('TITLE=%s',filenames{i}); end

    z = specs{i,4};
    if z>0
        pm = (specs{i,3}-1.0072763+z*1.0072763)/z; 
        fprintf(fid,'BEGIN IONS\r\nPEPMASS=%.5f\r\nCHARGE=+%.0f\r\n',pm,z);
    else 
        pm = specs{i,3}; 
        fprintf(fid,'BEGIN IONS\r\nPEPMASS=%.5f\r\n',pm);
    end;
    if ~isempty(scan) fprintf(fid,'%s\n',scan); end;
    if ~isempty(title) fprintf(fid,'%s\n',title); end;
    if ~isempty(seq) fprintf(fid,'%s\n',seq); end;
    if ~isempty(specs{i,2}) fprintf(fid,'%.5f %.5f\r\n',double(specs{i,2}(:,1:2)')); end;
    fprintf(fid,'END IONS\r\n\r\n');
end;

fclose(fid);
