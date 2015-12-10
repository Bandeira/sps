function plotCoverage(protLen, peptides, pixelInterval, colors, blists, matchedRegions)
% function plotCoverage(protLen, peptides, pixelInterval, colors, blists, matchedRegions)
%
% Example: plotCoverage(178,[10 20; 15 25],10,[],[]);
%
%  colors - column vector of color intensities [0:1] for every entry in peptides
%  blists - lists of b-ion peaks present per peptide
%  matchedRegions - highlight regions of the protein recovered (i.e. grouped/sequenced). List of start/end positions as nx2 vector.
%
peptidesIdx = find(peptides(:,1)>0);   peptides = peptides(peptidesIdx,:);   
if ~isempty(colors) colors = colors(peptidesIdx,:); end;   if ~isempty(blists) blists=blists(peptidesIdx,:); end;

numPeptides = size(peptides,1);

plot([1:protLen], pixelInterval*(numPeptides+1)*ones(1,protLen),'-');
hold on;

[foo idxS] = sortrows(peptides);   peptides = peptides(idxS,:);   
if ~isempty(colors) colors = colors(idxS,:); end;   if ~isempty(blists) blists=blists(idxS,:); end;

if size(colors,2)==1
    idxPos = find(colors>=.5); idxNeg = find(colors<.5);
	colorsPN = ones(size(colors,1),3);   % [colors zeros(size(colors,1),1) 1-colors];
	colorsPN(idxNeg,3) = 2*(.5-colors(idxNeg));   colorsPN(idxPos,1) = 2*(colors(idxPos)-.5);
else colorsPN=colors; end;

h = pixelInterval*(numPeptides+1);
if ~isempty(matchedRegions) matchedRegions = matchedRegions(find(matchedRegions(:,2)>0),:); end;
for r=1:size(matchedRegions,1)
%     idx = find(peptides(:,1)>=matchedRegions(r,1) & peptides(:,2)<=matchedRegions(r,2));
%     if isempty(idx) idx = find((peptides(:,2)>matchedRegions(r,1) & peptides(:,1)<=matchedRegions(r,1)) | (peptides(:,1)<=matchedRegions(r,2) & peptides(:,2)>=matchedRegions(r,2))); end;
%     if isempty(idx) idx = [1 size(peptides,1)]; end;  % if no spectra intersect the matched region then mark it from top to bottom
%     rectangle('Position',[matchedRegions(r,1) h-max(idx)*pixelInterval matchedRegions(r,2)-matchedRegions(r,1) (max(idx)-min(idx)+1)*pixelInterval],'FaceColor',[1 .8 0],'LineWidth',0.1);
% %     else
% %         rectangle('Position',[matchedRegions(r,1) .95*h max(matchedRegions(r,2)-matchedRegions(r,1),1) .05*h+2*pixelInterval],'FaceColor',[.9 0 0],'LineWidth',0.1);
% %     end
    rectangle('Position',[matchedRegions(r,1) 0 max(matchedRegions(r,2)-matchedRegions(r,1),1) h],'FaceColor',[.8 .8 .8],'LineWidth',0.1);
end

for p=1:numPeptides
    yPos = pixelInterval*(numPeptides-p+1);
if peptides(p,1)==420
    a=1;
end;
    if (~isempty(colors))
        plot([peptides(p,1):peptides(p,2)], yPos*ones(1,peptides(p,2)-peptides(p,1)+1),'-','Color',colorsPN(p,:));
    else
        plot([peptides(p,1):peptides(p,2)], yPos*ones(1,peptides(p,2)-peptides(p,1)+1),'-');
    end;
    if ~isempty(blists) & ~isempty(blists{p}) plot([peptides(p,1)+blists{p}(:,1)], yPos*ones(size(blists{p},1),1), '+r'); end;
end;
axis([0 protLen 0 numPeptides+2]);
hold off;
