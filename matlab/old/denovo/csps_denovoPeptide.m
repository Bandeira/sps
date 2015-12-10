function strDenovo = denovoPeptide(spec, szJumps, in_masses, in_letters, peakTol)
% function strDenovo = denovoPeptide(spec, szJumps, in_masses, in_letters, peakTol)
%
%  Converts every consecutive peak distance to an amino acid mass
%

global AAmasses AAletters;
if nargin<2 szJumps=1; end;
if nargin<4 szJumps=1; in_masses=AAmasses; in_letters=AAletters; end;
if nargin<5 szJumps=1; in_masses=AAmasses; in_letters=AAletters; peakTol=0.5; end;

if szJumps==1 [jumps, letters] = csps_getjumps(in_masses, in_letters, szJumps);
else [jumps, letters] = csps_csps_getjumpsAll(in_masses, in_letters, szJumps, .01); end;

strDenovo = '';
numPeaks = size(spec,1);   if numPeaks==0 return; end;
if spec(1,1)<=0.0001 masses = spec(2:numPeaks,1)-spec(1:numPeaks-1,1);
else masses = spec(:,1)-[0; spec(1:numPeaks-1,1)]; end;
for m=1:length(masses)
    idx = find(abs(jumps-masses(m))<=peakTol);
    if ~isempty(idx)
        if szJumps==1
            if length(idx)>1 aaStr=sprintf('[%s',letters{idx(1)}); else aaStr=sprintf('%s',letters{idx(1)}); end;
            for i=2:length(idx) t=sprintf('%s,%s',aaStr,letters{idx(i)}); aaStr=t; end;
            if length(idx)>1 t=sprintf('%s]',aaStr); aaStr=t; end;
        else aaStr=mergeStrings(letters(idx)); end;
        strDenovo = sprintf('%s%s',strDenovo,aaStr);
    else strDenovo = sprintf('%s[%.2f]',strDenovo,masses(m)); end;
end;

function s = mergeStrings(strs)
% Merges a set of lettersSep entries (as returned by getJumpsAll)

single = {};   for i=1:size(strs,1) single = [single; strs{i}]; end;
sizes = zeros(size(single,1),1);   for i=1:size(single,1) sizes(i)=size(single{i},2); end;
single = single(find(sizes==min(sizes)));

if size(single,1)==1 s=single{1}; else
    s=sprintf('[%s',single{1});
    for i=2:size(single,1) str=sprintf('%s,%s',s,single{i}); s=str; end;
    str=sprintf('%s]',s); s=str;
end;
