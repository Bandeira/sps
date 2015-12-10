function sn_save_abinfo(abinfo, filename)
% function sn_save_abinfo(abinfo, filename)
%
%  Saves information about the ABruijn components created by masab.
%
%  abinfo - 2-col cell, one line per component:
%         Col 1: 2-col int array with [specIndex flipState (0-no flip, 1-flip)]
%         Col 2: cell vector, 1 line per ABruijn vertex. Each ABruijn vertex is a 2-col vector of [specIndex peakIndex].
% 
% File format is:
%    0 (zero) [int]
%    Major version number [2 byte unsigned short]
%    Minor version number [2 byte unsigned short]
%    numUsedSpectra [unsigned int]
%    numUsedSpectra-by-2 unsigned int array of [used spectrum index,component index]
%    numUsedSpectra-by-1 1-byte array of [specFlipped = 0/1]
%    numComponents [unsigned int]
%    numComponents-by-1 unsigned short array of number of ABruijn vertices per component
%    numABVertices-by-1 unsigned short array of number of spectrum peaks per ABruijn vertex
%    totNumSpecPeaks-by-1 unsigned int array of spectrum index per peak (per ABruijn vertex)
%    totNumSpecPeaks-by-1 float array of peak masses per spectrum peak

fid = fopen(filename,'w'); if fid<=0 fprintf(1,'Error opening %s!\n',filename); return; end;
szUsed = 0;   fwrite(fid,szUsed,'uint');
szUsed = 1;   fwrite(fid,szUsed,'ushort');  % Version 1.1
szUsed = 1;   fwrite(fid,szUsed,'ushort');

used = [];   % spectrum index (col.1), flipped/not-flipped (col.2), component index (col.3)
numComps = size(abinfo,1);   
cInfo = zeros(numComps,1);  % #ABruijn vertices per component
vInfo = [];  % #peaks per ABruijn vertex
peaks = [];  % all ABruijn/spectrum-peak memberships: spectrum index (col.1), peak mass (col.2)
for c=1:numComps
    used=[used; abinfo{c,1} repmat(c,size(abinfo{c,1},1),1)];
    cInfo(c) = size(abinfo{c,2},1);
    peakCounts = zeros(cInfo(c),1);   cPeaks=[];
    for v=1:cInfo(c)
        peakCounts(v) = size(abinfo{c,2}{v},1);
        if peakCounts(v)>0
            cPeaks = [cPeaks; [abinfo{c,2}{v}(:,1)-1 abinfo{c,2}{v}(:,2)] ];
        end
    end
    vInfo = [vInfo; peakCounts];
    peaks = [peaks; cPeaks];
end
used = sortrows(used);   szUsed = size(used,1);
fwrite(fid,szUsed,'uint');
fwrite(fid,used(:,[1 3])'-1,'uint');  % Indices should be saved as zero-based
fwrite(fid,used(:,2),'char');

fwrite(fid,numComps,'uint');
fwrite(fid,cInfo,'ushort');
fwrite(fid,vInfo,'ushort');
fwrite(fid,peaks(:,1),'uint');
fwrite(fid,peaks(:,2),'float');
