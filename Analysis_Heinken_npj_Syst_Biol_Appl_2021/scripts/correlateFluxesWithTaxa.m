% This script extracts the correlations between taxa and computed secretion 
% fluxes for the personalized microbiome models generated for the 
% PLEASE/COMBO cohort.
% Almut Heinken, 03/21

% Compute Spearman correlations between net secretion fluxes and taxa
analysedFluxPath = [pwd filesep 'NetFluxes' filesep 'netProductionFluxes.csv'];
[FluxCorrelations, PValues, TaxonomyInfo] = correlateFluxWithTaxonAbundance(abunFilePath, analysedFluxPath, 'AGORA_infoFile.xlsx', 'Spearman');

% Export results for plots
levels=fieldnames(FluxCorrelations);
for i=1:length(levels)
    data = FluxCorrelations.(levels{i});
    % extract only strong correlations
    % first rows
    cnt=1;
    inds=[];
    for j=2:size(data,1)
        vals=cell2mat(data(j,2:end));
        if ~any(abs(vals) >0.75)
            inds(cnt,1)=j;
            cnt=cnt+1;
        end
    end
    data(inds,:)=[];
    % then columns
    cnt=1;
    inds=[];
    for j=2:size(data,2)
        vals=cell2mat(data(2:end,j));
        if ~any(abs(vals) >0.75)
            inds(cnt,1)=j;
            cnt=cnt+1;
        end
    end
    data(:,inds)=[];
    
    data(1,2:end)=strrep(data(1,2:end),',',' ');
    
    cell2csv([corrPath filesep 'Correlations_',levels{i},'.csv'],data);
    if ~strcmp(levels{i},'Phylum')
    % extract taxonomy info
    createTaxonomy = TaxonomyInfo.(levels{i});
    [C,IA] = setdiff(createTaxonomy(:,1),data(2:end,1),'stable');
    createTaxonomy(IA(2:end),:)=[];
    cell2csv([corrPath filesep 'Taxonomy_Correlations_',levels{i},'.csv'],createTaxonomy);
    end
end
