function seqs = load_denovo_seqs(filenames_prefix)
% function seqs = load_denovo_seqs(filenames_prefix)
%
%  Loads a set of denovo reconstructions constructed with cpplib's denovo
%

seqs=[];
allSeqs = load_pklbin(strcat(filenames_prefix,'.pklbin'), 1);     if isempty(allSeqs) return; end;
index = load_binArray(strcat(filenames_prefix,'.index'), 'int');  if isempty(index) fprintf(1,'(warning) File not found: %s.index\n',filenames_prefix); seqs=allSeqs; return; end;

seqs = cell(size(index,1),1);   curCount=0;
for i=1:size(index,1)
    seqs{i} = allSeqs(curCount+1:curCount+index(i,2),:);
    curCount = curCount+index(i,2);
end
