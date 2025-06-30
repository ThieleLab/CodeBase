%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Master script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AUTHOR: Tim Hensen,
% INSTITUTE: University of Galway, School of medicine

% DESCRIPTION:
% This script contains the complete computational analysis described in 
% "Metabolic modelling links gut microbiota to metabolic markers of
% Parkinson's disease". All input materials can be freely downloaded from 
% the Harvard Dataverse: https://doi.org/10.7910/DVN/PUTXL9. Before running
% this code, please inputs in a folder "input" under the project directory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Investigate cohort metadata %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, the cohort metadata will be described.
clear;clc;
%initCobraToolbox
% Set paths for analysis
paths = struct;
paths.root = what('parkinson_recreated').path; % IMPORTANT: DEFINE THE PROJECT DIRECTORY HERE
paths.inputs = fullfile(paths.root,'inputs');
paths.outputs = fullfile(paths.root,'outputs'); % Set directory to save analysis results

paths.metadataPath = fullfile(paths.inputs,'metadata.xlsx'); % Set path to metadata

% Force delete all current outputs ... 
if isfolder(paths.outputs)
    rmdir(paths.outputs,'s');
end

% and create a new clean folder for the
% outputs.
mkdir(paths.outputs)


% First, we will process the metadata file. Specifically, the column name
% of sample_ID will be converted to ID and male/female indications will be
% converted from m/f to male/female.
metadata = processMetadata(paths.metadataPath);

% The metadata file contains 724 samples, including 490 PD patients and 234
% controls. However, not all samples include all required metadata. In the
% next function, metadata variables of interest are defined and samples
% with one or more missing information in these metadata variables are
% removed. We also removed samples that were obtained using the swap
% method. The metadata variables of interest were obtained from Wallen et
% al. (https://doi.org/10.1038/s41467-022-34667-x).
[prunedMetada, sampleNumbers] = pruneParkinsonMetadataSamples(metadata);

% Add processed and pruned metadata to paths variable
paths.metadataPath = fullfile(paths.outputs,'pruned_processsed_metadata.xlsx');
writetable(prunedMetada,paths.metadataPath);
writetable(sampleNumbers,fullfile(paths.outputs,'metadata_removed_sample_statistics.txt'))

% After pruning, 654 samples were used for all further analysis, including
% 435 PD patients and 219 controls. Next, we will produce a table
% summarising all sameple metadata.

% Create a summary statistics table for the metadata
stats = makeMetadataTable(prunedMetada);
writetable(stats,fullfile(paths.outputs,'metadataSummaryStats.xlsx'),'WriteRowNames',true)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Metagenomics read count mapping %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Metagenomics read count mapping was performed against the APOLLO+AGORA2
% resources using the MARS pipeline. Note that mapping was performed on the
% unfiltered dataset with 724 samples. MARS also performed descriptions
% ot the microbiomes before and after mapping. Paths to the MARS mapping
% results can be found here:
%
% Number of microbial species and total read counts before and after
% mapping and the mapping coverages.
MARS_outputs = struct;
MARS_outputs.mappingStats = fullfile(paths.root, 'resultMARS','metrics','Species','summ_stats_Species.csv');
% Phylum abundances before mapping
MARS_outputs.PhylumPreMapping = fullfile(paths.root,'resultMARS','metrics','Species','postMapping_present_abundanceMetrics_Phylum.csv');
% Phylum abundances after mapping
MARS_outputs.PhylumPreMapping = fullfile(paths.root,'resultMARS','metrics','Species','preMapping_abundanceMetrics_Phylum.csv');
% Species abundances before mapping
MARS_outputs.SpeciesPreMapping = fullfile(paths.root,'resultMARS','metrics','Species','preMapping_abundanceMetrics_Species.csv');
% Species abundances post mappnig
MARS_outputs.SpeciesPreMappingAbundance = fullfile(paths.root,'resultMARS','metrics','Species','postMapping_absent_abundanceMetrics_Species.csv');


% Prepare gut microbial species results from Wallen 2022
paths.parkinsonMicrobes = fullfile(paths.inputs,'PD_associatedTaxa_subsetPresentSpeciesMARS.csv'); % Taxonomic annotations of microbes

% Process gut microbiome data
paths.ProcessedParkinsonMicrobes = fullfile(paths.outputs,'wallen2022_speciesRes.csv');
speciesData = prepareWallenSpeciesResults(paths.parkinsonMicrobes, fullfile(paths.outputs,'wallen2022_speciesRes.csv'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Describe investigated metabolites                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The list of PD-associated metabolites of interest was curated from a list
% of 190 metabolites with replicated PD associations, of which 127 were
% present in the WBMs. These associations were found by Luo et al. (2024),
% https://doi.org/10.1038/s41531-024-00732-z.
% Next, we identified which metabolites were also present in the blood
% compartment of the WBMs, and which metabolites were present in the gut
% microbiota lumen and diet.

% Clean up workspace
clearvars -except paths; clc;

% Set inputs
paths.models = fullfile(paths.inputs,'models');
paths.mWbms = fullfile(paths.models,'HMmodels');
paths.pdMets = fullfile(paths.inputs,'luo_2024_replicated_PD_mets.xlsx');
paths.metOntology = fullfile(paths.inputs,'metaboliteChemicalAnnotation.xlsx');
    

% Check metabolites of interest for presence in the blood, microbiota, and
% diet compartments.
[metabolitesCheck, metaboliteSummaryStats] = checkMetabolitePresenceInWBMs(paths);

% Add metabolite location info path to paths variable
paths.metWbmLocations = fullfile(paths.outputs,'metaboliteWbmLocationInfo.csv');

% Save tables metabolite location info and summary statistics
writetable(metabolitesCheck, paths.metWbmLocations)
writetable(metaboliteSummaryStats, fullfile(paths.outputs,'metaboliteWbmLocationInfo_summary.csv'),'WriteRowNames',true)


% 116 metabolites were present in the WBM blood compartments.
% Then, we described the biochemical annotations of these metabolites and
% extracted which metabolic subsystems these metabolites were involved in.
[biomchemClass,capturedSubsystems] = biochemClassifications(paths.metOntology);

% Save biochemical ontology summary statistics to excel file with a
% sheet for each ontology level
arrayfun(@(x) writetable(biomchemClass.(x), ...
    fullfile(paths.outputs,'metaboliteOntologySummary.xlsx'), 'Sheet', x), ...
    string(fieldnames(biomchemClass)), 'UniformOutput', false);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Processing of FBA solution                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% After metagenomic mapping, microbiome community models were generated
% using mgPipe.m, after which host-microbiome WBMs were created using
% createBatchMWBM.m (PERSEPHONE). Harvey/Harvetta 1.04c were used for
% creating these WBMs. The average european diet was given to the models.
% Then, analyseWBMsPD.m was run to investigate maximal metabolite productions
% in the WBM blood compartments for the 116 metabolites. Here, we will
% perform custom steps for processing the FBA results that do not exist in
% PERSEPHONE.

% Clean up workspace
clearvars -except paths; clc;

% Set new folder path to outputs from flux processing
paths.fluxes = fullfile(paths.outputs,'fluxes');

% Create flux processing output folder
if ~isfolder(paths.fluxes); mkdir(paths.fluxes); end

% Remove all FBA results that were not in the metadata after sample
% pruning.
fluxPath = pruneFBAResults(fullfile(paths.inputs, 'resultFlux'), paths.fluxes, paths.metadataPath);

% Rename the old SpeciesBIO field to taxonNames 
OK = ensurePersephoneCompatibility(fluxPath);

% Next, the parameters for flux processing are set. Metabolites for which a
% solution could be found in 5% of samples or less are removed.
paramFluxProcessing.rxnRemovalCutoff = {'fraction',  0.05};
paramFluxProcessing.fluxMicrobeCorrelationMetric = 'spearman_rho';

% Process flux results with the following Persephone function
analyseWBMsolPD(fluxPath,paramFluxProcessing, paths.fluxes);  % Function can take long time

% Set path to the processed fluxese
paths.fluxPath = fullfile(paths.fluxes,'processed_fluxes.csv');

% After analysing the FBA results, we check for metabolites that are linearly
% dependent with an upstream metabolite in the same pathway. I will
% investigate these linear dependencies here:
[prunedFluxes, linearDepMetaboliteTable] = pruneLinearlyDependentResults( ...
    paths.fluxPath, ...
    paths.metWbmLocations, ...
    paths.fluxes);

% The linearly dependent fluxes will still be investigated separately to
% ensure the most honest multiple testing correction after statistical
% analysis.

% Add paths to the paths variable as inputs for further processing of the
% FBA results.
paths.fluxLimDir = fullfile(paths.fluxes,'biomass_shadow_prices'); % biomass prices for finding flux-associated microbes
paths.mappedMicrobePath = fullfile(paths.fluxes,'WBM_relative_abundances.csv'); % Microbial relative abundances

% Fix relative abundance and shadow price files by combining the results of
% Collinsella masiliensis and Enorma masiliensis, which are the same
% species.
updated = fixMassiliensisFiles(paths.mappedMicrobePath, paths.fluxLimDir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     Differential flux analysis                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% After having prepared the flux data, we will associated the predicted
% fluxes against PD status, while controlling for relevant confounders.

% Clean up workspace
clearvars -except paths; clc; close all;

% Set output path to figures
paths.figures = fullfile(paths.outputs,'figures');

% Create output folder for figures
if ~isfolder(paths.figures); mkdir(paths.figures); end

% Next, we investigate outliers in the flux results
close all; [sampleImportanceTable, explVar] = sampleOutliersInFluxes(paths.fluxPath, paths.metadataPath, paths.figures);

% We will now remove the top 3 samples in the sampleImportanceTable
paths.prunedMetadataPath = pruneFluxOutliersFromMetadata(sampleImportanceTable,paths.metadataPath, 3);

% Define the confounders to control for
confounders = {...
    'Age_at_collection',...
    'Sex',...
    'Do_you_drink_alcohol',...
    'Laxatives',...
    'Probiotic',...
    'Antihistamines',...
    'Depression_anxiety_mood_med',...
    'Pain_med',...
    'Sleep_aid',...
    'total_sequences',...
    'Antibiotics_current'};

% After processing the flux results, we associate the fluxes against PD
% status while controlling for the defined confounding factors

% Set name of output file
paths.parkinsonFluxes = fullfile(paths.outputs,'fluxRegressions.csv');

[regressionResults,preparedInputTable,preparedMetadata, ~, regressions] = performParkinsonAnalysis( ...
    paths.fluxPath, ...
    paths.metadataPath, ...
    confounders, ...
    paths.parkinsonFluxes);

% Get reactions of interest

% Visualise regression results. We will use an FDR
% threshold of 0.108 to also include myristic acid, which is an interesting
% metabolites
FDRthreshold = 0.108;
close all; createViolinsPD(regressionResults,preparedInputTable,preparedMetadata, paths.figures, FDRthreshold);

% Update the paths variable with new fields
paths.rxnsOfInterest = regressionResults.Flux.Reaction(regressionResults.Flux.FDR < FDRthreshold );

% Investigate the influence of age and sex on the regression results
[fullModRes, paths.modResPath] = moderationAnalysisPD(preparedInputTable, preparedMetadata, confounders, paths.rxnsOfInterest, paths.outputs);

% Investigate if lifestyle variables should be added. 

[lifeStyleInfluence,filePath] = analyseLifeStyleFactorsPD(paths, confounders, regressions, paths.outputs);
% No lifestyle variables will be added as confounders to the regressions
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Identifying flux-associated microbial species               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% After obtaining the flux results, we identify which microbial species
% potentially influenced the flux results using shadow price values of 
% pan microbe biomass metabolites.
% We used a cutoff of microbial species that influenced the
% flux result in at least 5% of the samples.

% Clean up workspace
clearvars -except paths; clc; close all;

% Set new folder path to outputs from flux processing
paths.microbeToFlux = fullfile(paths.outputs,'microbeToFlux');

% Create flux processing output folder
if ~isfolder(paths.microbeToFlux); mkdir(paths.microbeToFlux); end

% Set paths
paths.mContributionDir = fullfile(paths.fluxes,'metabolic_influence_scores');

% Find flux-associated microbial species
[fmAssociations, fmLinksSummary, paths.fmLinksPath, ~] = findFluxLimiters(paths.mContributionDir, paths.rxnsOfInterest, paths.microbeToFlux);

% Find the microbial species with the highest average flux contributions
bootSamp = 50000;
[bootMeanTable, ~] = fluxMicrobeSensitivityAnalysis(paths.fmLinksPath, paths.mContributionDir, paths.microbeToFlux, bootSamp);
%%
% Filter on the most important microbial species and add the 
% path to the most important species to the paths variable.
cumulativeFraction = 0.95;
[microbialContributors, bootMeanCellFull, paths.filteredfmLinksSummaryPath,paths.mContributionTablePath] = filterMostImportantMicrobes(bootMeanTable, cumulativeFraction, paths.microbeToFlux);

% Find the flux-associated microbial species again for the filtered list
% and find summary statistics
[fmAssociationsFiltered, fmLinksSummaryFiltered, ~, filePathStats] = findFluxLimiters(paths.mContributionDir, paths.rxnsOfInterest, paths.microbeToFlux, microbialContributors);

% % To create the upset plot, you can run this function from R
% % C:\Users\mspg\Documents\wbm_modelingcode\src\ParkinsonAnalysis\functions\plotFigures.R
% 
% % Calculate jaccard distance. OPTIONAL. NOT DESCRIBED IN MANUSCRIPT
% [reorderedTable, meanJDist] = analyseSubsetSimilarities('', false,intersectionTable);
% 
% % Load the total mapped microbial speciespaths.mappedMicrobePath  OPTIONAL. NOT DESCRIBED IN MANUSCRIPT
% [phylumOverrep, annotatedSpecies] = phylumFluxAssocOverRepAnalysis(influencerSubset, paths.mappedMicrobePath, paths.parkinsonMicrobes, paths.outputs);
% close all % Close figures

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Flux-microbe correlation analysis                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace
clearvars -except paths; clc; close all;

% After having identified the flux-associated microbial species, we will
% correlate their relative abundances with the associated flux predictions.

% Correlate fluxes with microbe relative abundances. Then, we will prune
% the correlation table on microbial species that are in the top X
% strongest correlations for any of the metabolites.
[corrTable, corrTableTopMicrobes, filePaths] = correlateFluxWithRa(paths.fluxPath, paths.filteredfmLinksSummaryPath, paths.mappedMicrobePath, paths.outputs);
paths.fluxRaCorrPath = filePaths.fluxCorrPath;
paths.prunedFluxCorrPath = filePaths.prunedFluxCorrPath;

% Now perform an enrichment analysis
if 0 % The enrichment analysis is de-prioritized, but will be added later. 
    resolution = 'Phylum'; % Taxonomic level of interest
    saveDir = paths.outputs; % Folder of saved results
    % Perform enrichment analysis on the top correlated microbial species
    % taxonEnrichmentAnalysis
    [enrichmentRes, taxonomyTable, pathToResults] = taxonEnrichmentAnalysis(paths, resolution, corrCutoff);
end

% After finding the top microbial correlates for the flux predictions, we
% will find clusters of microbial species that better predict the fluxes.

% Set function parameters
top=5;
onlyCorrTopMicrobes = false;
cutoff = 0.05;

tic 
% Perform cluster correlations
[bestCorrPerSizeList, minbestCorrPerSizeList, rxnsOfInterest] = clusterFluxRaCorr(paths, top,onlyCorrTopMicrobes);
toc % ~100 seconds for top 4

% Add the correlation coefficients of the summed relative abundances for
% the sets of flux-associated microbes
bestCorrPerSizeList = addCorrForSummedFluxLinkedTaxa(corrTable,bestCorrPerSizeList, rxnsOfInterest);
minbestCorrPerSizeList = addCorrForSummedFluxLinkedTaxa(corrTable,minbestCorrPerSizeList, rxnsOfInterest);

% Concatenate correlations in a single table
topCorrelations = vertcat(bestCorrPerSizeList{:},minbestCorrPerSizeList{:});
topCorrelations = removevars(topCorrelations,'clusterIndex');

% Save results to table
paths.topClusterCorr = fullfile(paths.microbeToFlux, 'top_microbialClusterCorrelates.csv');
writetable(topCorrelations,paths.topClusterCorr)

% Generate adjacency matrix for visualisation for the metabolite-microbe
% assocations.
[adjMatrix,filePath] = extractTopFluxMicrobeLinks(paths);


%% Supplementary tables

% Get the supplementary tables 
% Clean up workspace
clearvars -except paths; clc;

% Table 1: Mapped and unmapped taxa + mapping coverages

% Load the present and absent species and concatinate the results. 
% Add a one for all the mapped species and a zero for all the unmapped
% ones. Then, add taxonomic and phylogenetic information.
clc
suplFilePath = fullfile(paths.outputs,'Supplementary_tables.xlsx');
if isfile(suplFilePath); delete(suplFilePath); end % For script testing

% Preset variable with table information


%%%% MICROBIOTA AND MAPPING STATISTICS ->

% Set inputs for species-level summary statistics
paths.species = fullfile(paths.inputs,'resultMARS', 'metrics','Species'); % Set directory with MARS results
fileNames = {'preMapping_abundanceMetrics_Species.csv',... % Define the files to load
    'postMapping_present_abundanceMetrics_Species.csv',...
    'postMapping_absent_abundanceMetrics_Species.csv'};

% Generate supplementary table
mergedTaxa = collectAbundanceStats(paths, fileNames);


% Set inputs for phylum-level summary statistics
microbiotaPath = fullfile(paths.inputs,'resultMARS','normalized','normalized_species.csv'); 
taxonomyPath = fullfile(paths.inputs,'resultMARS','preprocessedInput_afterRenaming.csv');
microbiotaWbmPath = fullfile(paths.fluxes,'WBM_relative_abundances.csv');

% Collect summary statistics on phylum-level mapping effects
summaryMerged = collectPhylumStats(microbiotaPath, taxonomyPath, microbiotaWbmPath);

% Save summary statistics to table
description = cell(2,1); 
description{1} = 'Mean relative abundances of mapped and unmapped microbial phyla and species'; % Header
description{2} = ['Left table: Microbial species mapping status and mean (SD) relative abundances across the cohort. ',...
    'Right table: Microbial phylum-level mean (SD) relative abundances across the cohort for the total, mapped, and unmapped microbial species.']; % Details
folder = paths.outputs;
suplTable = mergedTaxa;
writeSupplement(mergedTaxa, description, paths.outputs, summaryMerged)


%%%% Microbiome-WBM content statistics ->

% Get model summary statistics
disp('Start model content analysis')
if 1 % This part can take 168 seconds with 18 parallel workers on an intel core i9-10890xe
    tic
    [modelStats, summaryStats] = getPdMicrobiomeWBMstats(paths.mWbms, feature('numCores')); % 
    toc 
    
    
    % Save model content information to table
    description = cell(2,1); 
    description{1} = 'Summary statistics of microbiome-WBM model content'; % Header
    description{2} = ''; % Details
    writeSupplement(summaryStats, description, paths.outputs)
end
%%
%%%% DIET ->

% Find model diet and save to supplementary tables
dietTable = collectDietInfo(paths.mWbms);
% Save dietary information to table
description = cell(2,1); 
description{1} = 'Diet makeup of host-microbiome WBMs'; % Header
description{2} = 'Dietary metabolites taken from the average European diet (DOI: http://dx.doi.org/10.1093/nar/gky992)'; % Details
writeSupplement(dietTable, description, paths.outputs)


%%%% METABOLITE INFO ->

% Combine metabolite ontology information with metabolite names and
% presence in WBM compartments
suplMetaboliteTable = collectMetInfo(paths.metWbmLocations,paths.metOntology);
% Save metabolite information to table
description = cell(2,1); 
description{1} = 'Analysed metabolites with replicated metabolomic PD associations'; % Header
description{2} = 'The PD associated metabolites, VMH mappings, and biochemical annotations were taken from Luo et al. (2024), https://doi.org/10.1038/s41531-024-00732-z'; % Details
writeSupplement(suplMetaboliteTable, description, paths.outputs)

%%%% REGRESSION RESULTS ->

preparedfluxRegressionTable = prepareFluxRegressionResults(paths.parkinsonFluxes); % Flux regression results
preparedfluxModerationRegressionTable = prepareFluxModerationRegressionResults(paths.modResPath); % Flux regression results for moderation of age and sex

% Save dietary information to table
description = cell(2,1); 
description{1} = 'Logistic regression outcomes of predicted fluxes against Parkinson diagnosis'; % Header
description{2} = ['The logistic regressions were performed on log2-transformed and z-scaled blood fluxes. ', ...
    'The regression log odds represent the estimated change in log odd probability of PD with an increase in predicted log2 flux of one standard deviation.', ...
    'Positive regression coefficients indicate positive associations between predicted fluxes and PD, while negative regression coefficients indicate negative correlations with the fluxes.', ...
    'The regression log odds for the age-flux interaction represent the estimated change in the regression coefficient for flux and PD status with an increase in age of 1 log2 z-scaled unit of age in years', ...
    'The regression log odds for the sex-flux interaction represent the estimated change in the regression coefficient for flux and PD status with a change from male to female samples.'];% details
writeSupplement(preparedfluxRegressionTable, description, paths.outputs, preparedfluxModerationRegressionTable)


%%%% METABOLITE FLUX CONTRIBUTIONS ->

[mContributions,filteredfmLinksSummary] = collectFluxContributionInfo(paths.mContributionTablePath,paths.filteredfmLinksSummaryPath);

% Save microbial contribution summary to table
description = cell(2,1); 
description{1} = 'Microbial contributions to selected metabolites'; % Header
description{2} = ['Left table: Bootstrapped mean average and 95% confidence interval (CI) of microbial flux contributions in mmol/day/person. 50,000 boot samples were used to obtain the summary statistics. The last column indicates if a flux-microbe pair was included in further analyses. The smallest group of microbial species that contributed to 95% of the total flux contributions was analysed for each metabolite.',...
    'Right table: Summary of analysed metabolite-microbe pairs.']; % Details
writeSupplement(mContributions, description, paths.outputs, filteredfmLinksSummary)



%%%% FLUX-MICROBE CORRELATIONS ->

% Load flux-microbe correlations
fluxMicrobeCorr = readtable(paths.prunedFluxCorrPath,'VariableNamingRule','preserve','ReadRowNames',false);
fluxMicrobeCorr = renamevars(fluxMicrobeCorr,'Row','Microbial species');

% Save PD microbe information to table
description = cell(2,1); 
description{1} = 'Spearman correlations of analysed flux-microbe pairs'; % Header
description{2} = 'Spearman correlations were performed on the microbial relative abundances against predicted blood fluxes in mmol/day/person for each of the shown metabolites. Empty cells indicate microbial species that were excluded from this analysis due to having potential microbial flux contributions below the cutoff.';
writeSupplement(fluxMicrobeCorr, description, paths.outputs)


%%%% Microbial PD associations ->

% Load processed microbe file
microbePdInfo = readtable(paths.ProcessedParkinsonMicrobes, 'VariableNamingRule','preserve');

% Save PD microbe information to table
description = cell(2,1); 
description{1} = 'Gut microbial associations with PD diagnosis'; % Header
description{2} = 'Microbial species in the host-microbiome WBMs and their previously found associations with PD diagnosis. This table is adapted from results by Wallen et al. (2022), https://doi.org/10.1038/s41467-022-34667-x'; % Details
writeSupplement(microbePdInfo, description, paths.outputs)



%%% MICROBIAL CLUSTER CORRELATES -> 

% Load table
topClusterCorr = readtable(paths.topClusterCorr);

% Save cluster correlation results to supplementary table
description = cell(2,1); 
description{1} = 'Top flux-correlating microbial combinations'; % Header
description{2} = 'Spearman correlation coeficients of microbial combinations that best correlated with their associated flux predictions';
writeSupplement(topClusterCorr, description, paths.outputs)


% Move supplementary tables to sharepoint

sourcePath = fullfile(paths.outputs,'Supplementary_tables.xlsx');
targetPath = 'C:\Users\mspg\National University of Ireland, Galway\Group_MSP - Parkinson\Supplementary_tables.xlsx';

if isfile(targetPath)
    delete(targetPath)
end

copyfile(sourcePath,targetPath)
