
abunFilePath = [rootDir filesep 'normalizedCoverage.csv'];
fluxPath = [rootDir filesep 'Modeling_CRC' filesep 'Plots' filesep 'DrugFluxes_CRC_Microbiomes.csv'];
infoFilePath = 'AGORA2_infoFile.xlsx';

%% calculate correlations between fluxes and taxon abundance-Japanese diet
[FluxCorrelations, PValues] = correlateFluxWithTaxonAbundance(abunFilePath, fluxPath, infoFilePath, 'Spearman');
cell2csv([rootDir filesep 'Modeling_CRC' filesep 'Plots' filesep 'JD_Species.csv'],FluxCorrelations.('Species'));
close all

%% Get taxonomy for R plot
createTaxonomy=readInputTableForPipeline('AGORA2_infoFile.xlsx');
% remove columns not needed
findTaxonLevel=find(strcmp(createTaxonomy(1,:),'Species'));
createTaxonomy(:,1:findTaxonLevel-1)=[];
% reduce to only unique entries
[C,IA,IC] = unique(createTaxonomy(:,1),'stable');
createTaxonomy=createTaxonomy(IA,1:size(createTaxonomy,2));
% find entries not in data
[C,IA] = setdiff(createTaxonomy(:,1),FluxCorrelations.('Species')(2:end,1),'stable');
createTaxonomy(IA(2:end),:)=[];
createTaxonomy(:,10:end)=[];
cell2csv([rootDir filesep 'Modeling_CRC' filesep 'Plots' filesep 'Taxonomy_JD_Species.csv'],createTaxonomy);

