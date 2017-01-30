
function [nonConflicAssoc]= nonConflictingAssociations(jstar)

global T psi

nonConflicAssoc= [];

for j= 1:psi
    if T(j,end) > T(jstar,end)
        if ~any( ( T(jstar,1:end-1) - T(j,1:end-1) ) .* T(j,1:end-1))
            nonConflicAssoc= [nonConflicAssoc, j];
        end
    end
end







