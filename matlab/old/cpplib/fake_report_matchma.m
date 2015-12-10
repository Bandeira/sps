function [report, reportScores] = fake_report_matchma(specs, specType, specsName, mods, modMasses, addModMasses, peakTol, bypercs)
%
%  Uses the peptide annotations in specs to create a report in matchma format
%

if isempty(bypercs)
	stage82(specs, specType, peakTol, mods, modMasses, addModMasses, 0, 0, specsName, 0);
	eval(sprintf('load %s_stage82 %s_bypercs; bypercs = %s_bypercs; clear %s_bypercs;',specsName,specsName,specsName,specsName));
end;

numSpecs = size(specs,1);   report = cell(numSpecs,1);   if size(specs,2)==5 idxPep=5; else idxPep=7; end;
for i=1:numSpecs
    report{i} = cell(1,13);
    if bypercs(i,1)>=bypercs(i,3) report{i}{1,5}='b'; report{i}{1,9}=bypercs(i,5); 
    else report{i}{1,5}='y'; report{i}{1,9}=bypercs(i,7); end;
    report{i}{1,8} = specs{i,idxPep};
end

reportScores = [max(bypercs(:,[5 7])')' max(bypercs(:,[1 3])')'];
