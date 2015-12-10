function report_signatures(name, mincos, topk, specplotPath)
% function report_signatures(name, mincos, topk)
%
%  Creates cytoscape and cluster_info files to report results for spectral
%  networks analysis of metabolites MS/MS spectra
%

specs = sn_load_pklbin('spectra/specs_ms.pklbin');
txt=repmat({''},size(specs,1),1);   txt2 = txt;
for i=1:length(txt)
    pm = specs{i,3};
    if specs{i,4}>0
        pm = (pm + specs{i,4}-1) / specs{i,4};  % Convert back to precursor mass
    end;
    txt{i}=sprintf('s%d\t%.3f',i,pm); 
    if length(specplotPath)>0
        txt2{i}=sprintf('%s/specplot --pklbin spectra/specs_ms.pklbin --spectrum %d --outdir png --outfile spectrum_%d.png --title "Cluster spectrum %d, precursor mass %.3f"',specplotPath,i,i,i,pm);
    end;
end;
sn_save_lines(sprintf('%s_cytoscape_parent_masses.txt',name),txt,1);
if length(specplotPath)>0
    sn_save_lines('genpng.sh',['mkdir png'; txt2],1);
end

clstinfo = read_csv('ClusterInfo.csv',';',4);
clstinfo = [clstinfo repmat({''},size(clstinfo,1),1)];
txt = repmat({''},size(clstinfo,1),1);
for i=1:size(clstinfo,1)
    j = str2num(clstinfo{i,1})+1;
    txt{i} = sprintf('%d\t%s\t%s\t%s\t%.3f',j,clstinfo{i,2},clstinfo{i,3},clstinfo{i,4},specs{j,3});
end
sn_save_lines(sprintf('%s_cluster_info.txt',name),txt(:,1),1);

pairs = sn_load_aligns('aligns/pairs_raw.bin',6);
fprintf(1,'%7d pairs loaded\n',size(pairs,1));
pairs = pairs(find(pairs(:,5)>=mincos),:);
fprintf(1,'%7d pairs with cos >= %.3f\n',size(pairs,1),mincos);
pairs_topk = filterpairs_topk(pairs,topk,2);
fprintf(1,'%7d pairs in top %d for both vertices\n',size(pairs_topk,1),topk);
f=fopen(sprintf('%s_cytoscape_top%d_both.txt',name,topk),'w');
for i=1:size(pairs_topk,1) fprintf(f,'s%d s%d %.4f\n',pairs_topk(i,1),pairs_topk(i,2),pairs_topk(i,5)); end;
fclose(f);
