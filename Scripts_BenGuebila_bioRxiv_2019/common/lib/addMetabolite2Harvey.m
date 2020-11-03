function harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met,stochvec)
% Lumps reaction flux constraints through adding dummy metabolite
% Input    : harvey    - COBRA model structure
%            oldNumVec - index of reactions fluxes lumped together
%            organ     - organ name
%            met       - metabolite name
%            OPTIONAL 
%            stochvec - vector of +1 or -1 for reaction directionality correspondance
% Output   : New model structure

%fix reactions through adding dummy metabolite
if nargin<5
    if isequal(met,'glc_D[c]')
        stochvec= harvey.A(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_' met])),harvey.mets)),oldNumVec);
    elseif isequal(met,'g6p[c]')
        stochvec= harvey.A(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_' met])),harvey.mets)),oldNumVec);
    end
end
vec2Add = vecBuilder(oldNumVec,harvey,stochvec);
%add dummy metabolite
harvey.S = [harvey.S;vec2Add];
harvey.S = sparse(harvey.S);
%add dummy met in b vector
harvey.b = [harvey.b;0];
harvey.mets = [harvey.mets;'GIM_dummy'];
%add met names, formula, charge
harvey.metNames = [harvey.metNames;'GIM_dummy'];
harvey.metFormulas = [harvey.metFormulas;'GIM_dummy'];
%harvey.metCharge = [harvey.metCharge;0];

end