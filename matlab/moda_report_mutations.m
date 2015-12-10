function moda_unq = moda_report_mutations(fn_msgfdb, fn_moda, fn_moda_unq, ppmTol)
% function moda_unq = moda_report_mutations(fn_msgfdb, fn_moda, fn_moda_unq, ppmTol)
%
% Extracts putative mutated peptide IDs from MODa IDs using the following
% filters:
%  1 - remove IDs for scans identified by MS-GFDB
%  2 - remove IDs whose modification can be explained as a common PTM
%  3 - enforce ppm precursor mass tolerance
%
% Filtered results are returned in moda_unq and saved in tab-separated
%  format to fn_moda_unq
%

global AAmasses AAletters AtomicMasses;
massH2O = [0 2 0 1 0] * AtomicMasses;

% Load raw results
msgfdb = read_csv(fn_msgfdb,'\t', 19);  msgfdb = msgfdb([2:size(msgfdb,1)],:);   szM = size(msgfdb,1);
moda =  read_csv(fn_moda,'\t', 16);     
szA = size(moda,1)-1;
moda_head = [moda(1,:) {'BaseAA'} {'MutMass'} {'MutAA'}];
moda = [moda([2:szA+1],:) cell(szA,3)]; % Extra 3 cells: BaseAA, MutMass, MutAA

% Get MODa-only results
all = cell(szM+szA,6); % concat filename:scan, M/A, original index, 
                       %  MODa unique=0/1 or MSGFDB=2, filename (col.5),
                       %  scan (col.6)
for i=1:szM
    all(i,:) = [{strcat(msgfdb{i,1},msgfdb{i,3},'.M')} {'M'} {i} {2} msgfdb(i,[1 3])];
end
for i=1:szA
    all(szM+i,:) = [{strcat(moda{i,1},moda{i,3},'.A')} {'A'} {i} {0} moda(i,[1 3])];
end
[~,idxS] = sort(all(:,1));
all = all(idxS,:);
for i=1:szA+szM-1
    if all{i,2}=='A'
        t = strcat(all{i,1}(1:length(all{i,1})-2),'.M');
        if strcmp(t, all{i+1,1})==0 % Not identified by MS-GFDB, keep it
            all{i,4} = 1;
        end
    end
end
moda_unq = moda( [ all{ [all{:,4}]'==1 ,3} ], :);
save moda_unq1 moda_unq;
% load moda_unq1;
sz = size(moda_unq,1);
fprintf(1,'Got %d IDs from MODa with no MS-GFDB ID\n', sz);


%
% Ignore IDs with common mods: mass-delta/valid-sites/site-loc-tolerance
%
commonMods = [ [{'+16'} {'MW'}]; [{'-17'} {'.Q'}]; [{'-18'} {'*'}]; [{'+1'} {'*'}]; ; [{'+57'} {'.C'}] ];
siteLocTol = 2;
numMods = size(commonMods,1);
% sort mods by decreasing string length
commonMods = [commonMods repmat({'0'},numMods,1)];
for i=1:numMods
    commonMods{i,3} = length( commonMods{i,1} );
end
[foo,idxS]=sort([commonMods{:,3}]');   commonMods = commonMods(idxS(numMods:-1:1),:);
% find/clear unmodified peptides or peptides with common mods
for i=1:sz
    peptide = moda_unq{i,10}(3:length(moda_unq{i,10})-2);
    modPos = find( peptide =='+' | peptide =='-' );
    if isempty(modPos)
        moda_unq{i,10} = '';  % <<-- delete unmodified peptides
    else
        modStr = peptide(find(peptide<'A' | peptide>'Z'));
        modIdx = find( strcmp( modStr, commonMods(:,1) ) );
        if ~isempty(modIdx)
            if commonMods{modIdx,2}=='*'
                moda_unq{i,10} = '';  % <<-- delete common mod, not site-specific
            else
                for siteIdx=1:length(commonMods{modIdx,2})
                    site = commonMods{modIdx,2}(siteIdx);
                    if site=='.'
                        if modPos<=siteLocTol
                            moda_unq{i,10} = '';  % <<-- delete common mod, N-term
                            break;
                        end
                    else
                        siteOffset = find(peptide==site) - modPos;
                        if ~isempty(siteOffset) && ~isempty( find( siteOffset>=-siteLocTol-1 & siteOffset<=siteLocTol+commonMods{modIdx,3} ) ) % commonMods{modIdx,3} = mod string length
                            moda_unq{i,10} = '';  % <<-- delete common mod, site within range
                            break;
                        end
                    end
                end  % Looping over sites
            end
        end % Done processing mods
    end
end
moda_unq = moda_unq( find( strcmp( moda_unq(:,10),'' ) == 0)  ,:);
save moda_unq2 moda_unq;
% load moda_unq2;
sz = size(moda_unq,1);
fprintf(1,'Got %d IDs from MODa after removing unmodified IDs or common mods\n', sz);

% Select only mutation offsets: mass-delta/valid-sites/site-loc-tolerance
intAAMasses = round(AAmasses);
for i=1:sz
    peptide = moda_unq{i,10}(3:length(moda_unq{i,10})-2);
    modPos = find( peptide =='+' | peptide =='-' );
	if modPos>1
        modStr = peptide(find(peptide<'A' | peptide>'Z'));
        modSite = peptide(modPos-1);
        aaMass = intAAMasses( AAletters == modSite );
        modMass = str2double(modStr);
        
        idx = find( intAAMasses == (aaMass+modMass) );
        if isempty(idx)
            moda_unq{i,10} = '';  % <<-- delete non-mutation peptides
        else
            pepMass = sum(sn_getmasses( peptide(find(peptide>='A' & peptide<'Z')) ,'',[],0)) + massH2O;
            expMass = str2double(moda_unq{i,4});
            match = 0;
            for j=1:length(idx)
                thMass = pepMass - AAmasses( AAletters == modSite ) + AAmasses( idx(j) );
                if 1000000*abs(expMass-thMass)/thMass <= ppmTol
                    moda_unq{i,17} = modSite;
                    moda_unq{i,18} = sprintf('%.4f',expMass - pepMass);
                    moda_unq{i,19} = AAletters(idx(j));
                    match = 1;
                    break;
                end
            end
            if match == 0 
                moda_unq{i,10} = '';  % <<-- delete non-mutation peptides
            end
        end
    end
end
moda_unq = moda_unq( find( strcmp( moda_unq(:,10),'' ) == 0)  ,:);
save moda_unq3 moda_unq;
% load moda_unq3;
sz = size(moda_unq,1);
fprintf(1,'Got %d mutation IDs from MODa with ppm error <= %d\n', sz, ppmTol);

% Find mutations covered by 2+ unique peptides --> same prot, same mod, peptide pos overlap

% save moda_unq as tsv
save_tsv(fn_moda_unq,[moda_head; moda_unq]);
