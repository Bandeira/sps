function sn_grid_prep(numSpecs, numJobs, exeDir, commonParams)
%
% Warning: assumes that each qsub has access to 2 CPUs
%

numSlices=numJobs*2;   pairs = [(numSpecs-1):-1:0]';   pairs = cumsum(pairs);   indices = ones(numSlices,1);  % i-th pos: first spectrum index for the i-th job
loadsize = round(pairs(size(pairs,1))/numSlices);      jobSize = zeros(numSlices,1);
for i=2:numSlices
    idx = find(pairs>=loadsize)+1;   idx=idx(find(idx<=length(pairs)));   % Have each job compute a little over its load quota to avoid concentrating leftovers on the last job
    if ~isempty(idx) indices(i) = max(min(idx),indices(i-1)+1);
    else
        numSlices = i-1; indices = indices(1:numSlices);
        fprintf(1,'Warning: number of jobs reduced to %d\n',numJobs);
        break;
    end;
    jobSize(i-1) = pairs(indices(i)-1);
    pairs = pairs-pairs(indices(i)-1);
    loadsize = round(pairs(size(pairs,1))/(numSlices-i+1));  % Each job takes a load quota a little over its share and thus changes the load size for the remaining jobs by a little bit
end
indices = [indices; numSpecs+1];  % Last entry in indices should be first index to ignore

fn_means = cell(2*numJobs,1);   fn_vars = cell(2*numJobs,1);   fn_aligns = cell(2*numJobs,1);   fn_ratios = cell(2*numJobs,1);
fidMain = fopen(sprintf('run_jobs.sh'),'w');
for job=1:numJobs
    jobStr = sprintf('%2d',job);   jobStr(find(jobStr==' '))='0';   % replace all spaces with zeros to enforce correct filename sorting (i.e. 01<10)
    fidJob = fopen(sprintf('job%s.sh',jobStr),'w');
    fprintf(fidJob,'#!/bin/bash\n#\n#$ -cwd\n#$ -j y\n#$ -S /bin/bash\n#\n\n');

    for slc=1:2
        sliceIdx = 2*(job-1)+slc;   if sliceIdx > numSlices break; end;
        sliceStr = strcat(jobStr,sprintf('%d',slc));
        fid = fopen(sprintf('job%s.params',sliceStr),'w');
        fprintf(fid,'%d\n',6+size(commonParams,1));
        fprintf(fid,'IDX_START=%d\nIDX_END=%d\n',indices(sliceIdx)-1,indices(sliceIdx+1)-2);  % -2 to avoid overlaps between end/starts of consecutive jobs
        fprintf(fid,'OUTPUT_ALIGNS=aligns_%s.bin\nOUTPUT_MEANS=means_%s.bin\nOUTPUT_VARIANCE=vars_%s.bin\nOUTPUT_RATIOS=ratios_%s.bin\n',sliceStr,sliceStr,sliceStr,sliceStr);
        for p=1:size(commonParams,1) fprintf(fid,'%s\n',commonParams{p}); end;
        fclose(fid);
        
        fprintf(fidJob,'%s/filterpairs job%s.params &\n',exeDir,sliceStr);
        fn_means{sliceIdx} = sprintf('means_%s.bin',sliceStr);    fn_vars{sliceIdx} = sprintf('vars_%s.bin',sliceStr);
        fn_aligns{sliceIdx} = sprintf('aligns_%s.bin',sliceStr);  fn_ratios{sliceIdx} = sprintf('ratios_%s.bin',sliceStr);
    end
    fprintf(fidJob,'wait\n');    fclose(fidJob);
    fprintf(fidMain,'qsub -l h_vmem=200M job%s.sh\n',jobStr);
end
fclose(fidMain);

sn_save_lines('fn_means.txt', fn_means, 1);    sn_save_lines('fn_vars.txt', fn_vars, 1);
sn_save_lines('fn_aligns.txt', fn_aligns, 1);  sn_save_lines('fn_ratios.txt', fn_ratios, 1);
