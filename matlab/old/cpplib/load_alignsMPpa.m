function specsMP = load_alignsMPpa(specs, aligns, filename)
%
%  Loads the sets of matched peaks on each partial-overlap alignment
%    (as output by batch_pa)
%
%  specsMP(:,1) - sets of indices grouping matched-peaks-specs from each spectrum pair;
%                   has size(specsMP,1)/4 elements, one per pair in aligns
%  specsMP(i,2:5) - matched peaks (col.2), parent mass (col.3), unused (col.4), peptide (col.5)
%

specsMP = load_pkl(filename);   numSpecs=size(specsMP,1);   numAligns=size(aligns,1);   idxPep=size(specs,2);
if abs(numSpecs/4 - numAligns)>0.0001 specsMP=[]; fprintf(1,'ERROR: Number of spectra in %s does not match number of pairs in aligns!\n',filename); return; end;

pairs = [aligns(:,1:2) aligns(:,1:2)];   specsIdx = reshape(pairs',numSpecs,1);   specsMP(:,5) = specs(specsIdx,5);
sets = reshape([1:numSpecs]',4,numAligns)';   for i=1:numAligns specsMP{i,1}=sets(i,:); end;
