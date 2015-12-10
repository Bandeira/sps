function [prmList1,prmList2] = sepbpm(prms1, mhMass1, prms2, mhMass2, shift, tolerance)
% function [prmList1 prmList2] = sepbpm(prms1, mhMass1, prms2, mhMass2, shift, tolerance)
%
%  Separate the prms1 and prms2 into two sets of prms such that at most one prm per twin pair matches something
%    in the other spectrum. Uses maximum biparte matching - guarantess that no twin pairs are used but not a 
%    minimum 56 Da distance between PRMs.
%

[pairs1, scores1] = getPrmPairs(prms1, mhMass1, tolerance);   sz1 = size(pairs1,1);   scores1 = mean(scores1')';
pairs1=[pairs1 [1:sz1]'];   pairs1 = [[pairs1(:,[1 3]); pairs1(sz1:-1:1,[2 3])] [scores1; scores1(sz1:-1:1)]];

[pairs2, scores2] = getPrmPairs(prms2, mhMass2, tolerance);   sz2 = size(pairs2,1);   scores2 = mean(scores2')';
pairs2=[pairs2 [1:sz2]'];   pairs2 = [[pairs2(:,[1 3]); pairs2(sz2:-1:1,[2 3])] [scores2; scores2(sz2:-1:1)]];
pairs2(:,2) = pairs2(:,2)+sz1;  % Make the vertice numbers different 

pairs1(:,1) = pairs1(:,1)/10;   pairs2(:,1) = pairs2(:,1)/10+shift;

edges = [];
for i=1:size(pairs1,1)
    idx = find(abs(pairs2(:,1)-pairs1(i,1))<=2*tolerance);
    for j=1:size(idx,1) edges = [edges; pairs1(i,2) pairs2(idx(j),2) pairs1(i,3)+pairs2(idx(j),3)]; end;

%     idx = find(abs(pairs2(:,1)-shift-pairs1(i,1))<=2*tolerance);
%     for j=1:size(idx,1) edges = [edges; pairs1(i,2) pairs2(idx(j),2) pairs1(i,3)+pairs2(idx(j),3)]; end;
end

match = MaxMatch(edges);   

edges(:,2) = edges(:,2)-sz1;  % Recover the original prm pair numbers for spectrum 2
prmList1 = [];   prmList2 = [];
for i=1:size(match,1)
    p1 = edges(match(i),1);   p1 = [p1 size(pairs1,1)+1-p1];
    p2 = edges(match(i),2);   p2 = [p2 size(pairs2,1)+1-p2];
    
    for j=1:2
        for k=1:2
            if abs(pairs1(p1(j),1) - pairs2(p2(k),1)) <= 2*tolerance
                prmList1 = [prmList1; pairs1(p1(j),[1 3])];   
                prmList2 = [prmList2; pairs2(p2(k),1)-shift pairs2(p2(k),3)];
            end
        end
    end
%    res = [res; pairs1(p1,1) pairs2(p2,1) pairs1(size(pairs1,1)+1-p1,1) pairs2(size(pairs2,1)+1-p2,1)];
end

[foo idx] = sort(prmList1(:,1));   prmList1 = prmList1(idx,:);
[foo idx] = sort(prmList2(:,1));   prmList2 = prmList2(idx,:);

% TODO: impose sparse subset + generate same type of consensus for y-matches (denovoMA takes 2 prmLists as parameters)