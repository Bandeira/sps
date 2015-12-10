function specSet = load_pklbin(filename, convertToDouble)
% function specSet = load_pklbin(filename, convertToDouble)
%
%  Load a pklbin file as created by cpplib.  convertToDouble converts all spectral data to double (as usual)
%
%  specSet has 5 column format with index (col.1), spectra (col.2), parent masses (col.3) and parent charges (col.5)
%

if nargin<2 convertToDouble=1; end;
fid=fopen(filename,'r');  if fid<=0 fprintf(1,'Error opening file %s!\n',filename); specSet=[]; return; end;

numSpecs = fread(fid,1,'int32');   specSet = cell(numSpecs,5);
numPeaks = fread(fid,numSpecs,'int16=>double');
for i=1:numSpecs
    if convertToDouble data = fread(fid,2*numPeaks(i)+2,'float32=>double');
    else data = single(fread(fid,2*numPeaks(i)+2,'float32')); end;
    data2 = reshape(data,2,numPeaks(i)+1)';
    specSet{i,2} = data2(2:numPeaks(i)+1,:);
    specSet{i,3} = double(data2(1,1));
    specSet{i,4} = round(double(data2(1,2)));
%     specSet{i,5} = '';
end

fclose(fid);
