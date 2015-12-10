function [binIndex, fileInfo] = binPM(filenames, filetype, resolution, tolerance, maxMass, topBins, baseFN)
% function [binIndex, fileInfo] = binPM(filenames, filetype, resolution, tolerance, maxMass, topBins, baseFN)
%
%  Reads spectrum precursor masses pre-extracted (grep) from MS/MS spectrum files, bins them
%   for the given resolution and find the topBins (allowing for non-overlapping tolerance).
%
%  Each of the topBins is output to a separate text file baseFN_<bin>.txt with the following format:
%    <num lines>
%    <spectra filename>:<spec_index0>:...:<spec_indexN>
%
%  binIndex - Index of which files went into which bins. binIndex(i,:) - integer bin mass (col.1), ordered set of pairs: input spectrum file index (pos.1), spectrum index in input file (pos.2)
%  fileInfo(i,:) - file name (col.1), ordered list of scan info per spectrum (col.2), vector of spectrum precursor masses (col.3)
% 
%  Notes:
%   - Output spectrum indices are one-based
%   - Output <spectra filenames> always have .pklbin extensions
%   - filenames{i} is expected to be of the format <name>.filetype.txt
%

if ~strcmp(filetype,'ms2') fprintf(1,'Currently only supports ms2 formats\n'); return; end;

% Compute the bins (h) while keeping track of where each spectrum falls (index)
maxMass = ceil(maxMass/resolution);   h = zeros(maxMass,1);   index = cell(maxMass,1);
numFiles = size(filenames,1);   fileInfo = cell(numFiles,3);
if exist('binPM_data.mat') load binPM_data; else
    for fIdx=1:numFiles
        fprintf(1,'Processing file %d of %d...',fIdx,numFiles); tic;
        lines = load_lines(filenames{fIdx});   if isempty(lines) continue; end;
        
        % Extract precursor masses from ms2 files
        scansIdx = [1:3:size(lines,1)]';    scans = lines(scansIdx);
        linesIdx = [2:3:size(lines,1)]';    lines = lines(linesIdx);   numLines = size(lines,1);
        masses = zeros(numLines,1);         
        for i=1:numLines 
            mass = str2num(lines{i}(1:min(find(lines{i}==' '))-1));    masses(i) = mass;
            mass = round( mass/resolution );   if mass>maxMass mass=maxMass; end;
            h(mass) = h(mass)+1;   index{mass} = [index{mass}; [single([fIdx i])] ];
        end;
        filenames{fIdx} = filenames{fIdx}(1:length(filenames{fIdx})-4-length(filetype));
        
        fileInfo(fIdx,:) = [ filenames(fIdx) {scans} {masses} ];
        fprintf(1,'done in %f seconds.\n',toc);
	end
    save binPM_data h index filenames fileInfo;
end;

% Compute the tolearnce-smoothed hTol version of h
intTol = round(tolerance/resolution);   tolRange = [-intTol:intTol];
if exist('hTol')~=1
	hTol = zeros(size(h));
	indices = [1+intTol:size(h,1)-intTol]';
	indicesH = repmat(indices,1,size(tolRange,2)) + repmat(tolRange,size(indices,1),1);
	hTol(indices) = sum( h(indicesH)' )';
	save binPM_data h hTol index;
end;

% Greedily select the topBins
% [foo, idxS] = sort(hTol);   idxS = idxS(size(hTol,1):-1:1);
hTol(size(hTol,1)-intTol:size(hTol,1)) = -1;  % The last bins contain spectra with a mixture of parent masses > maxMass and should be ignored
chosenBins = zeros(topBins,1);   chosenIdx=1;   numChosen=0;   hTol = [hTol; [1:intTol]'];
[maxVal, maxIdx] = max(hTol);
if exist('binPM_chosenBins.mat') load binPM_chosenBins; else
	while maxVal>0
        maxIdx = maxIdx(1); 
        chosenBins(chosenIdx) = maxIdx;   chosenIdx = chosenIdx+1;   numChosen=numChosen+hTol(maxIdx); if chosenIdx>topBins break; end;
        
        % Avoid choosing overlapping bins - selected bin centers should be at least 2*intTol+1 away, i,j st. abs(i-j)>=2*intTol+1
%         for i=maxIdx-intTol:maxIdx+intTol hTol(i+tolRange)=hTol(i+tolRange)-h(i); end;
        hTol(maxIdx-2*intTol:maxIdx+2*intTol)=-1;
	
        [maxVal, maxIdx] = max(hTol);
	end
	if chosenIdx<topBins chosenBins = chosenBins(1:chosenIdx-1,:); end;
	fprintf(1,'Got %d bins with a total of %d spectra (%.1f average #spectra/bin)\n',size(chosenBins,1),numChosen,numChosen/size(chosenBins,1));
	save binPM_chosenBins chosenBins;
end

% Output the selection results - version.1, too slow - opens _all_ input files to extract only a few spectra
if isempty(baseFN) return; end;
% for i=1:size(chosenBins,1)
%     lclSpecIndex = zeros(sum(h(chosenBins(i)+tolRange)),2);   lclIdx=1;
%     for rIdx = chosenBins(i)-intTol:chosenBins(i)+intTol if h(rIdx)>0 lclSpecIndex(lclIdx:lclIdx+h(rIdx)-1,:) = double(index{rIdx});  lclIdx = lclIdx+h(rIdx); end; end;
%     lclSpecIndex = sortrows(lclSpecIndex);
%     
%     % Generate the output lines
%     filesIdx = unique(lclSpecIndex(:,1));
%     lclLines = [{sprintf('%d',length(filesIdx))}; cell(length(filesIdx),1)];
%     for fIdx = 1:length(filesIdx)
%         specsIdx = lclSpecIndex(find(lclSpecIndex(:,1)==filesIdx(fIdx)),2);
%         lclLines{1+fIdx} = sprintf('%spklbin:%d',filenames{filesIdx(fIdx)},length(specsIdx));
%         for sIdx=1:length(specsIdx) str=sprintf('%s:%d',lclLines{1+fIdx},specsIdx(sIdx)); lclLines{1+fIdx}=str; end;
%     end
%     save_lines(sprintf('%s_%d.txt',baseFN,chosenBins(i)), lclLines, 1);
% end

% Output the selection results - version.2, each input pklbin is open only once and its spectra are exported to the bins one at a time
%  output index file is of the format
%    <number of input files>
%    <input file name>:<indices list size>:<which bin_index takes spectrum 0>:...:<which bin_index takes spectrum n>
%
%  Notes:
%   - correspondence between bin_index and a bin file name is established with a companion <baseFN>_index.txt
%   - bin_index values are one-based, bin_index==0 means that the corresponding spectrum is not to be output
%
chosenBins = sort(chosenBins);   numBaseDigits = ceil(log10(length(chosenBins)));   numBins = size(chosenBins,1);
lines = [{sprintf('%d',length(chosenBins))}; cell(length(chosenBins),1)]; for lIdx=1:length(chosenBins) lines{lIdx+1}=sprintf('%s_%d.ms2',baseFN,chosenBins(lIdx)); end;
save_lines(sprintf('%s_filenames.txt',baseFN), lines, 1);

hIdx  = reshape(repmat(chosenBins',size(tolRange,2),1)+repmat(tolRange',1,numBins), size(chosenBins,1)*size(tolRange,2), 1);
hBins = reshape(repmat([1:length(chosenBins)],size(tolRange,2),1), numBins*size(tolRange,2), 1);
[hIdx,uIdx] = unique(hIdx);   hBins = hBins(uIdx);   if length(uIdx)~=numBins*size(tolRange,2) fprintf(1,'Warning: bins contained repeated elements!\n'); end;

% Create global index of which bin contains which spectra
specIndex = zeros(sum(h(hIdx)),3);   lclIdx=1;   % specIndex(i,:) - spectrum file index (col.1), spectrum index in file (col.2), bin index (col.3)
for rIdx = 1:length(hIdx) if h(hIdx(rIdx))>0 specIndex(lclIdx:lclIdx+h(hIdx(rIdx))-1,:) = [double(index{hIdx(rIdx)}) hBins(rIdx)*ones(h(hIdx(rIdx)),1)];  lclIdx = lclIdx+h(hIdx(rIdx)); end; end;
specIndex = sortrows(specIndex);     uInput = unique(specIndex(:,1));

binIndex = cell(numBins,2);  % binIndex(i,:) - integer bin mass (col.1), ordered set of pairs: input spectrum file index (pos.1), spectrum index in input file (pos.2)
binCount = zeros(numBins,1); % counts how many spectra were already added to a given bin
for bIdx=1:size(binIndex,1) binIndex{bIdx,1} = chosenBins(bIdx);   binIndex{bIdx,2} = zeros(size(find(specIndex(:,3)==bIdx),1),2); end;

% Create output index file and binIndex
lines = [{sprintf('%d',length(uInput))}; cell(length(uInput),1)];
for rIdx = 1:length(uInput)
    lclSpecIndex = specIndex(find(specIndex(:,1)==uInput(rIdx)),:);
    str = sprintf('%spklbin:%d',filenames{uInput(rIdx)}(1:length(filenames{uInput(rIdx)})-7),max(lclSpecIndex(:,2)));  prevIndex=0;   strPos = length(str)+1;
    lines{rIdx+1} = sprintf('%s%s',str,repmat(' ',1,(numBaseDigits+1)*max(lclSpecIndex(:,2))) );  % Additional spaces reserve memory for the whole line
    for sIdx=1:size(lclSpecIndex,1)
        while lclSpecIndex(sIdx,2)>prevIndex+1 lines{rIdx+1}(strPos:strPos+1)=':0'; strPos=strPos+2; prevIndex=prevIndex+1; end;  % Create zero entries for the spectra that are not to be output to any bin
        str=sprintf(':%d',lclSpecIndex(sIdx,3)); lines{rIdx+1}(strPos:strPos+length(str)-1)=str; strPos=strPos+length(str); prevIndex=lclSpecIndex(sIdx,2);
        
        bIdx = lclSpecIndex(sIdx,3);   binCount(bIdx) = binCount(bIdx)+1;
        binIndex{bIdx,2}(binCount(bIdx),:) = [uInput(rIdx) lclSpecIndex(sIdx,2)];
    end
    lines{rIdx+1} = lines{rIdx+1}(1:strPos-1);
end
save_lines(sprintf('%s_index.txt',baseFN), lines, 1);
