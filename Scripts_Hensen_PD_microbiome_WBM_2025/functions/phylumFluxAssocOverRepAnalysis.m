function [phylumOverrep, annotatedSpecies] = phylumFluxAssocOverRepAnalysis(influencerSubset, mappedMicrobePath, pdShift, saveDir)
% Phylum over-representation analysis of alll flux-associated microbiome
% subsets.


mappedSpecies = getSpeciesInWBMs(mappedMicrobePath);
% microbiome = readtable(mappedMicrobePath,'VariableNamingRule','preserve');
% mappedSpecies = table();
% mappedSpecies.Taxon = microbiome.Properties.VariableNames(2:end)';
% mappedSpecies.Taxon = replace(mappedSpecies.Taxon,'_',' ');
% mappedSpecies(end,:) = []; % Remove sum of taxa row

% Load table with mapped phylum information
taxaInfo = readtable(pdShift,'VariableNamingRule','preserve');
taxaInfo = renamevars(taxaInfo,'Species_AGORA2APOLLOtaxonomy','Taxon');

% Map the phylum info onto the mappedSpecies
taxaInfo = taxaInfo(:,{'Phylum','Taxon'});
annotatedSpecies = outerjoin(mappedSpecies, taxaInfo,'Type','left','Keys','Taxon','MergeKeys',true);

% Manually add the missing phylum information
annotatedSpecies.Phylum(matches(annotatedSpecies.Taxon,'Ruminococcus chamellensis')) = {'Firmicutes'};
annotatedSpecies.Phylum(matches(annotatedSpecies.Taxon,'Streptococcus anginosus')) = {'Firmicutes'};
annotatedSpecies.Phylum(matches(annotatedSpecies.Taxon,'Pseudomonas putida')) = {'Proteobacteria'};
annotatedSpecies.Phylum(matches(annotatedSpecies.Taxon,'Pseudomonas aeruginosa')) = {'Proteobacteria'};
annotatedSpecies.Phylum(matches(annotatedSpecies.Taxon,'Amedibacillus dolichus')) = {'Firmicutes'};
annotatedSpecies.Phylum(matches(annotatedSpecies.Taxon,'Comilactobacillus farciminis')) = {'Firmicutes'};

% Now add information to each microbial speces on its influence on the
% metabolite

% Define metabolites
metabolites = string(erase(influencerSubset(:,1),'.csv'));
% Map metabolic influencers on the total microbiome dataset
for i=1:length(metabolites)
    associatedMicrobes = erase(influencerSubset{i,2},'pan');
    associatedMicrobes = replace(associatedMicrobes,'_',' ');
    annotatedSpecies.(metabolites(i)) = matches(annotatedSpecies.Taxon,associatedMicrobes);
end

% Find all unique
phyla = unique(annotatedSpecies.Phylum);

% Preallocate overrepresentation table
phylumOverrep = cell(length(phyla),2);
phylumOverrep(:,1) = phyla;

for i=1:length(phyla)
annotatedSpecies.currPhylum = matches(annotatedSpecies.Phylum,phyla(i));
for j=1:length(metabolites)
    [conttbl,~,~,labels] = crosstab(annotatedSpecies.currPhylum,annotatedSpecies.(metabolites(j)));
    [~,p,stats] = fishertest(conttbl);
    phylumOverrep{i,j+1} = p;
end
end

phylumOverrep = [ [{'Phylum'} cellstr(metabolites')] ; phylumOverrep];


% Save results
writecell(phylumOverrep,[saveDir filesep 'phylumMicrobiomeSubsetEnrichment.csv'])
end