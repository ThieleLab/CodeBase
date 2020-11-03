function A = createDeltaMatchMatrixMod(set1,rxnsList)
%createDeltaMatchMatrix Create a flux difference constraint matrix for MOMA
%type calculations
% Marouen BEN GUEBILA 08/16 based on
% Markus Herrgard 1/4/07

nRxns1 = length(set1);

ind1 = rxnsList;

nCommon = length(rxnsList);

A = sparse(2*nCommon,nRxns1+2*nCommon);
for i = 1:nCommon
    A(i,ind1(i)) = -1;
    A(i,nRxns1+i) = 1;
    A(nCommon+i,ind1(i)) = 1;
    A(nCommon+i,nRxns1+nCommon+i) = 1;
end