
% This script contains the code for executing the genome-scale 
% reconstruction, microbiome model generation, simulation, and data
% analysis pipeline in APOLLO. Note that the generation of 
% genome-scale reconstructions and microbiome models, and execution of 
% simulations is extremely time-consuming (multiple weeks to months). Hence, the raw 
% simulation data will be provided to enable reproducing the analyses.

addpath('Scripts')

%% Refinement of draft genome-scale reconstructions and computation of model properties
% Draft genome-scale reconstructions created through KBase will be refined
% in the following steps (time scale: ~4-8 weeks)

% Pasolli genomes
runDemeter_Pasolli

% Almeida genomes
runDemeter_Almeida

%% Retrieval of taxonomic composition

plotTaxonomicComposition
getAPOLLOAGORA2Overlap

%% Analysis of model properties

% compute and plot various reconstruction statistics for Figure 2
getReconstructionStatistics
getMetaboliteStats
makeReconstructionSizeBoxplots
plotCorePanReactome

%% Creation and interrogation of microbiome models
% Personalized microbiome models will be created from relative strain-level 
% abundance data for 14,451 samples. For increased computational
% efficiency, the models are built in batches of 1,000 samples (time scale
% ~2 weeks).

% Body site microbiomes
prepareAbundanceDataBodySites
getModelsBodySites
startMgPipeBodySites

% Gut microbiomes
prepareAbundanceDataGut
getModelsGut
startMgPipeGut
collectMicrobiomeModels

% combine all data
combineDatasets

%% Definition of microbiome datasets and extraction of results
% The data for 11 defined datasets of microbiome samples will be extracted.
% For a subset of models, the production potential for some metabolites
% will also be computed (time scale: ~1 day).

analyseSubsets

%% Generation of summary tables on model properties
% From the computed model properties, summary tables and panels for Figure
% 2 will be generated. The taxonomic distribution will also be plotted for
% Figure 2.

%% Results from the machine learning classifier on strain level
% Extracts the results form the machine learning classifier on the
% strain-level reconstruction data. Creates a summary table and the input
% data for the panels in Figure 3.

extractRandomForestsResultsStrains
extractClassifyingFeatures

%% Results from the machine learning classifier on microbiome level
% Extracts the results from the machine learning classifier on the
% personalised microbiome modeling reconstruction data. Creates a summary 
% table and the input data for figures.

extractRandomForestsResultsMicrobiomes
extractClassifyingFeaturesMicrobiomes

%% Plotting of statistical analysis results for the 11 microbiome datasets
% Extracts the results from statistical analyses and creates a summary
% table and panels for Figure 5.

plotStatisticalAnalysisResults
extractSignificantData
plotSignificantSubsystems
extractStatisticalAnalysisResults
