function [V,E] = specc1p_rnd(szM)
%
%  Generate a matrix with the c1p property and with rows/cols permuted. 
%    Return a list of edges in aligns format (and a set of vertices in contigs format)
%

V = cell(szM,6);   for i=1:szM V{i,6}=i; end;
M = eye(szM);  for i=1:szM-1 M(i,i+1)=1; end;
P = randperm(szM);
M = M(:,P);   M = M(P,:);
% M = M(:,[1:szM/2 szM:-1:szM/2+1]);
M=M+M';  % make M symmetrical
[vi,vj]=find(M);
E = [vi vj ones(size(vi,1),2)];
