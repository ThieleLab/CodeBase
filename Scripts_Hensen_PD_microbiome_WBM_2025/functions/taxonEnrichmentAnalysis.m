function [enrichmentRes, taxonomyTable, filePath] = taxonEnrichmentAnalysis(paths, resolution, corrCutoff, saveDir)
%  INPUT:
% microbeAnnotationPath = paths.pdShift; % Path to microbe annotation data
% fluxCorrPath = filePaths.fluxCorrPath;
% corrCutoff = 0.1;

microbeAnnotationPath = paths.pdShift; % Path to microbial taxonomy info
fluxCorrPath = paths.fluxRaCorrPath; % Path to microbial correlations
saveDir = paths.saveDir; % Path to directory for saving the outputs

% Load and process microbe annotations
taxaInfo = readtable(microbeAnnotationPath,'VariableNamingRule','preserve');
taxaInfo.Species = replace(taxaInfo.Species,' ','_');
taxaInfo = taxaInfo(:,{'Phylum','Class','Order','Family','Genus','Species'});

% Load correlation results
corrInfo = readtable(fluxCorrPath,'VariableNamingRule','preserve','ReadRowNames',false);
corrInfo = renamevars(corrInfo,'Row','Species');

% Remove Flux-associated taxa

% Set all correlations under the threshold to zero
corrArray = corrInfo{:,2:end};
corrArray(abs(corrArray)<corrCutoff) = 0;
corrArray(corrArray~=0) = 1;
corrInfo{:,2:end} = corrArray;

% Merge taxa info with correlations
taxonomyTable = outerjoin(taxaInfo,corrInfo,'Type','right','Keys','Species','MergeKeys',true);

% Remove all non-taxa
taxonomyTable(matches(taxonomyTable.Species,{'Flux-associated taxa','Sum of taxa'}),:)=[];

% Manually add the missing phylum information
taxonomyTable.Phylum(matches(taxonomyTable.Species,'Ruminococcus_chamellensis')) = {'Firmicutes'};
taxonomyTable.Phylum(matches(taxonomyTable.Species,'Streptococcus_anginosus')) = {'Firmicutes'};
taxonomyTable.Phylum(matches(taxonomyTable.Species,'Pseudomonas_putida')) = {'Proteobacteria'};
taxonomyTable.Phylum(matches(taxonomyTable.Species,'Pseudomonas_aeruginosa')) = {'Proteobacteria'};
taxonomyTable.Phylum(matches(taxonomyTable.Species,'Amedibacillus_dolichus')) = {'Firmicutes'};
taxonomyTable.Phylum(matches(taxonomyTable.Species,'Comilactobacillus_farciminis')) = {'Firmicutes'};


% Get all unique taxa
taxa = table2cell(taxonomyTable(:,1));
uniqueTaxa = unique(taxa);

% Get all metabolites
uniqueMetabolites = corrInfo.Properties.VariableNames(2:end);

% Preallocate taxaEnrichmentRes table
numRows = length(uniqueMetabolites) * length(uniqueTaxa);
enrichmentRes = table('Size', [numRows, 8], ...
                'VariableTypes', {'string', 'string', 'double', 'double', ...
                                  'double', 'double', 'double', 'double'}, ...
                'VariableNames', {char(resolution), 'Metabolite', 'PValue', ...
                                  'taxonWithMetabolite', 'taxonWithoutMetabolite', ...
                                  'NonTaxonWithMetabolite', 'NonTaxonWithoutMetabolite','OR'});

% Counter for enrichmentRes table rows
rowIdx = 1;

% Perform Fisher's exact test for each genus and metabolite
for i = 1:length(uniqueTaxa)
    currTaxon = uniqueTaxa{i};
    taxonPresence = matches(taxa, currTaxon);
    
    for j = 1:length(uniqueMetabolites)
        currMetabolite = uniqueMetabolites{j};
        metabolitePresence = taxonomyTable.(currMetabolite);
        
        % Calculate contingency table values
        a = sum(taxonPresence & metabolitePresence);          % Metabolite with taxon
        b = sum(taxonPresence & ~metabolitePresence);         % No metabolite with taxon
        c = sum(~taxonPresence & metabolitePresence);         % Metabolite without taxon
        d = sum(~taxonPresence & ~metabolitePresence);        % No metabolite without taxon
        
        % Perform Fisher's exact test
        [~, pValue, stats] = fishertest([a, b; c, d], 'Tail', 'right');
        
        % Populate the enrichmentRes table
        enrichmentRes(rowIdx, :) = {currTaxon, currMetabolite, pValue, a, b, c, d, stats.OddsRatio};
        rowIdx = rowIdx + 1;
    end
end

% Add for each taxon if it associated more often with metabolite A than
% without metabolite A.
enrichmentRes.taxonPresence = enrichmentRes.taxonWithMetabolite ./ (enrichmentRes.taxonWithMetabolite+enrichmentRes.NonTaxonWithMetabolite);
enrichmentRes.taxonAbsence = enrichmentRes.taxonWithoutMetabolite ./ (enrichmentRes.taxonWithoutMetabolite + enrichmentRes.NonTaxonWithoutMetabolite);
enrichmentRes.Enrichment = enrichmentRes.taxonPresence ./ enrichmentRes.taxonAbsence;
% Add FDR adjustment
enrichmentRes.FDR = mafdr(enrichmentRes.PValue,'BHFDR',true);

% sort rows on p-value
enrichmentRes = sortrows(enrichmentRes,'PValue');

% Save enrichmentRes
filePath = [saveDir filesep append(char(resolution),'microbialEnrichmentsPD.csv')];
writetable(enrichmentRes,filePath)
end