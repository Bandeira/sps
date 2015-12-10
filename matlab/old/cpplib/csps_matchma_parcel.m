function parcels = matchma_parcel(str)
% function parcels = matchma_parcel(str)
%
%  Divides a string into {}, (), [] blocks and individual matched amino acids.
%
%  parcels - as many lines as there are parcels, each line is a parcel string (1), AA match start/end pos (2-3)
%

idxOpen = find(str=='{' | str=='(' | str=='[');
idxClose = find(str=='}' | str==')' | str==']');
numBlocks=length(idxOpen);   blocksLen = sum(idxClose-idxOpen+1);   numParcels = numBlocks+length(str)-blocksLen;

parcels = cell(numParcels,3);   curPos=1;   pIdx=1;   curAApos=1;
for bIdx=1:numBlocks
    if idxOpen(bIdx)>curPos
        numAA = idxOpen(bIdx)-curPos;
        for i=1:numAA
            parcels{pIdx,1}=str(curPos+i-1);   
            parcels{pIdx,2}=curAApos;   parcels{pIdx,3}=curAApos;   curAApos=curAApos+1;
            pIdx = pIdx+1;
        end
    end
    parcels{pIdx,1}=str(idxOpen(bIdx):idxClose(bIdx));
    if parcels{pIdx,1}(1)~='{' 
        parcels{pIdx,2}=curAApos;
        if parcels{pIdx,1}(1)=='(' parcels{pIdx,3}=curAApos+length(parcels{pIdx,1})-3; else
            idx = find(parcels{pIdx,1}==',');   

if length(idx)~=2 fprintf(1,'ERROR in matchma_parcel.m: Parcel string "%s" (from %s)does not have two commas (only %d) - exiting...\n',parcels{pIdx,1},str,length(idx)); return; end;

            txt = parcels{pIdx,1}(idx(1)+1:idx(2)-1);
            parcels{pIdx,3}=curAApos+length(txt)-1;
        end
        curAApos=parcels{pIdx,3}+1;
    end;
    pIdx = pIdx+1;   curPos = idxClose(bIdx)+1;
end

if isempty(idxOpen) | idxClose(numBlocks)<length(str)
    numAA = length(str)-curPos+1;
    for i=1:numAA
        parcels{pIdx}=str(curPos+i-1);
        parcels{pIdx,2}=curAApos;   parcels{pIdx,3}=curAApos;   curAApos=curAApos+1;
        pIdx = pIdx+1;
    end
end
