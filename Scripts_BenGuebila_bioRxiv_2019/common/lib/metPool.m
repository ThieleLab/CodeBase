function [harvey]=metPool(harvey,met)
    %Marouen BEN GUEBILA 08/2016
    %makes a pooled metabolite out of noni-intersecting metabolites
    a=harvey.S(cellfun(@(x) isequal(x,['RBC_' met]),harvey.mets),:) ;%
    harvey.S(cellfun(@(x) isequal(x,['RBC_' met]),harvey.mets),:) = 0;%delete
    c=harvey.S(cellfun(@(x)  isequal(x,['Platelet_' met]),harvey.mets),:) ;%
    harvey.S(cellfun(@(x)  isequal(x,['Platelet_' met]),harvey.mets),:) = 0;
    e=harvey.S(cellfun(@(x)  isequal(x,['Monocyte_' met]),harvey.mets),:) ;%
    harvey.S(cellfun(@(x)  isequal(x,['Monocyte_' met]),harvey.mets),:) = 0;
    g=harvey.S(cellfun(@(x)  isequal(x,['Nkcells_' met]),harvey.mets),:) ;%
    harvey.S(cellfun(@(x)  isequal(x,['Nkcells_' met]),harvey.mets),:) = 0;
    i=harvey.S(cellfun(@(x)  isequal(x,['CD4Tcells_' met]),harvey.mets),:) ;%
    harvey.S(cellfun(@(x)  isequal(x,['CD4Tcells_' met]),harvey.mets),:) = 0;
    k=harvey.S(cellfun(@(x)  isequal(x,['Bcells_' met]),harvey.mets),:) ;%
    harvey.S(cellfun(@(x)  isequal(x,['Bcells_' met]),harvey.mets),:) = 0;
    vec2Add = a+c+e+g+i+k;
    %add dummy metabolite
    harvey.S = [harvey.S;vec2Add];
    harvey.S = sparse(harvey.S);
    %add dummy met in b vector
    harvey.b = [harvey.b;0];
    harvey.mets = [harvey.mets;'GIM_dummy'];
    %add met names, formula, charge
    harvey.metNames = [harvey.metNames;'GIM_dummy'];
    harvey.metFormulas = [harvey.metFormulas;'GIM_dummy'];
end