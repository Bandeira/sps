function absShifts = getAbsoluteShifts(aligns, tolerance, printInfo)
% function absShifts = getAbsoluteShifts(aligns, tolerance, printInfo)
%
%  Finds absolute shifts in a directed connected graph where all shifts are coherent
%
%  absShifts - absolute shift from the start. Vector referenced by contig index
%

vertices = unique(aligns(:,[1,2]));
maxVertex = max(vertices);

absShifts = -inf*ones(maxVertex ,1);
counts    = zeros(size(absShifts));  % counts how many edges voted for a given shift (for averaging)

absShifts(vertices(1)) = 0;  % choose an arbitrary start
[absShifts, counts] = getAbsoluteShifts_aux(absShifts, counts, aligns, vertices(1), 0, tolerance, printInfo);

% uniqueness
shifts = absShifts(vertices);
[shiftsS, idxS] = sort(shifts);
[foo, idxBack] = sort(idxS);
cur = 1;
for i=2:size(shiftsS,1)
    if abs(max(shiftsS(cur))-shiftsS(i)) <= tolerance
        cur = [cur i];
    else
        shiftsS(cur) = sum(shiftsS(cur))/size(cur,2);
        cur = i;
    end
end
shiftsS(cur) = sum(shiftsS(cur))/size(cur,2);
shifts = round(10*(shiftsS(idxBack) - min(shiftsS)))/10;

absShifts(vertices) = shifts;


function [absShifts, counts] = getAbsoluteShifts_aux(absShifts, counts, aligns, startVertex, absVertexShift, tolerance, printInfo)
%
%  Note: originally a local function in mergeContigsSet.m

in = aligns(find(aligns(:,2)==startVertex),:);
out = aligns(find(aligns(:,1)==startVertex),:);

for i=1:size(in,1)
    incomingVertex = in(i,1);
    if absShifts(incomingVertex)==-inf
        absShifts(incomingVertex) = absVertexShift - in(i,3);
        counts(incomingVertex) = 1;
        [absShifts, counts] = getAbsoluteShifts_aux(absShifts, counts, aligns, incomingVertex, absShifts(incomingVertex), tolerance, printInfo);
    else
        absShifts(incomingVertex) = (counts(incomingVertex)*absShifts(incomingVertex) + absVertexShift - in(i,3)) / (counts(incomingVertex)+1);
        counts(incomingVertex) = counts(incomingVertex)+1;
        if printInfo & absShifts(incomingVertex) - (absVertexShift - in(i,3)) > tolerance;
            fprintf(1,'WARNING: Relative shifts do not match (delta = %.2f) !\n',absShifts(incomingVertex) - (absVertexShift - in(i,3)));
        end
    end
end
for i=1:size(out,1)
    nextVertex = out(i,2);
    if absShifts(nextVertex)==-inf
        absShifts(nextVertex) = absVertexShift + out(i,3);
        counts(nextVertex) = 1;
        [absShifts, counts] = getAbsoluteShifts_aux(absShifts, counts, aligns, nextVertex, absShifts(nextVertex), tolerance, printInfo);
    else
        absShifts(nextVertex) = (counts(nextVertex)*absShifts(nextVertex) + absVertexShift + out(i,3)) / (counts(nextVertex)+1);
        counts(nextVertex) = counts(nextVertex)+1;
        if printInfo & absShifts(nextVertex) - (absVertexShift + out(i,3)) > tolerance;
            fprintf(1,'WARNING: Relative shifts do not match (delta = %.2f) !\n',absShifts(nextVertex) - (absVertexShift + out(i,3)));
        end
    end
end
