function [specSet, info] = sn_load_pklbin(filename, convertToDouble)
% function [specSet, info] = sn_load_pklbin(filename, convertToDouble)

if nargin<2 convertToDouble=1; end;
fid=fopen(filename,'r');  if fid<=0 fprintf(1,'Error opening file %s!\n',filename); specSet=[]; return; end;

numSpecs = fread(fid,1,'int32');
if numSpecs>0
    specSet = load_pklbin_old(fid, numSpecs, convertToDouble);
    info = zeros(numSpecs,2);
else
    version = fread(fid,1,'char=>int');   subversion = fread(fid,1,'char=>int');
%     if version>2 || (version==2 && subversion>1)
%         fprintf(1,'Pklbin v%d.%d not supported.\n',version,subversion);
%     end
    if version==2 
        if subversion>1
            fprintf(1,'Pklbin v%d.%d not supported.\n',version,subversion);
            return;
        end
        [specSet, info] = aux_load_pklbin_2(fid, convertToDouble);
    end
    
    if version==3
        if subversion>1
            fprintf(1,'Pklbin v%d.%d not supported.\n',version,subversion);
            return;
        end
        [specSet, info] = aux_load_pklbin_3(fid, convertToDouble);
    end
    
end

fclose(fid);

function specSet = load_pklbin_old(fid, numSpecs, convertToDouble)
    specSet = cell(numSpecs,5);
    % fprintf(1,'DEBUG: got %d spectra\n',numSpecs);
    numPeaks = fread(fid,numSpecs,'int16=>double');
    for i=1:numSpecs
        if convertToDouble data = fread(fid,2*numPeaks(i)+2,'float32=>double');
        else data = fread(fid,2*numPeaks(i)+2,'float32'); end;
    % fprintf(1,'index %d: numPeaks(i) = %d, size(data)=[%d %d]\n',i,numPeaks(i),size(data,1),size(data,2));
        data2 = reshape(data,2,numPeaks(i)+1)';
        specSet{i,2} = data2(2:numPeaks(i)+1,:);
        specSet{i,3} = data2(1,1);
        specSet{i,4} = round(double(data2(1,2)));
        specSet{i,5} = '';
    end
 
function [specSet, info] = aux_load_pklbin_2(fid, convertToDouble)

    numSpecs = fread(fid,1,'int32');
    specSet = cell(numSpecs,5);
    scans = fread(fid,numSpecs,'uint32');
    msLevels = fread(fid,numSpecs,'int16=>double');
    fragType = fread(fid,numSpecs,'int16');
    parentMasses = fread(fid,numSpecs,'float=>double');
    parentCharge = fread(fid,numSpecs,'int16=>double');
    pmTols = fread(fid,numSpecs,'float');
    if subversion==0
        numPeaks = fread(fid,numSpecs,'int16=>double');
    else
        numPeaks = fread(fid,numSpecs,'uint=>double');
    end
    for i=1:numSpecs
        if convertToDouble data = fread(fid,3*numPeaks(i),'float32=>double');
        else data = fread(fid,3*numPeaks(i),'float32'); end;
        specSet{i,1} = scans(i);
        specSet{i,2} = reshape(data,3,numPeaks(i))';
        specSet{i,3} = parentMasses(i);
        specSet{i,4} = round(double(parentCharge(i)));
        specSet{i,5} = '';
    end
    info = [scans msLevels];

function [specSet, info] = aux_load_pklbin_3(fid, convertToDouble)

    numSpecs = fread(fid,1,'int32');
    specSet = cell(numSpecs,5);
    info = zeros(numSpecs,2);
    numStrings = fread(fid,1,'int32');
    for i=1:numStrings
        strSz = fread(fid,1,'int32');
        str = fread(fid,strSz,'char');
    end
    for i=1:numStrings
        value = fread(fid,1,'int16');
    end

    for specIdx=1:numSpecs
        scan = fread(fid,1,'uint32');
        msLevel = fread(fid,1,'int16=>double');
        info(specIdx,:) = [scan msLevel];
        fragmtType = fread(fid,1,'int16=>double');
        reversed = fread(fid,1,'char=>double');
        
        numStrings = fread(fid,1,'int32');
        for i=1:numStrings
            strSz = fread(fid,1,'int32');
            filename = fread(fid,strSz,'char');
        end;
        
        parentMass = fread(fid,1,'float=>double');
        parentMassTol = fread(fid,1,'float=>double');
        precCharge = fread(fid,1,'int16=>double');
        precMZ     = fread(fid,1,'double');

        numPeaks = fread(fid,1,'int32');
        masses   = fread(fid,numPeaks,'float=>double');
        tolerances   = fread(fid,numPeaks,'float=>double');
        intensities  = fread(fid,numPeaks,'float=>double');
        
        specSet{specIdx,1} = scan;
        specSet{specIdx,2} = [masses intensities tolerances];
        specSet{specIdx,3} = parentMass;
        specSet{specIdx,4} = round(precCharge);
        specSet{specIdx,5} = '';
    end
