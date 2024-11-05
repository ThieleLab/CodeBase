
% This script contains the code for executing the genome-scale 
% reconstruction, microbiome model generation, simulation, and data
% analysis pipeline in APOLLO. Note that the generation of 
% genome-scale reconstructions and microbiome models, and execution of 
% simulations is extremely time-consuming (multiple weeks to months). Hence, the raw 
% simulation data will be provided to enable reproducing the analyses.

addpath(genpath('code'))

%% MODEL BUILDING AND DATA GENERATION %%
% Run the following scripts to reproduce simulations (time-consuming)

%% Refinement of draft genome-scale reconstructions and computation of model properties
% Draft genome-scale reconstructions created through KBase will be refined
% in the following steps (time scale: ~8 weeks)

% download and extract draft reconstructions

% Pasolli genomes
runDemeter_Pasolli

% Almeida genomes
runDemeter_Almeida

%% Analysis of model properties

% compute various reconstruction statistics for Figure 2
getReconstructionStatistics
% compute model properties for machine learning analysis
computeModelProperties
combineAllStrainData

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
analyseMicrobiomes

%% DATA ANALYSIS %%
% Run the following scripts to analyse the raw simulation data

mkdir([rootDir filesep 'results'])

%% Retrieval of taxonomic distribution for Figure 2
plotTaxonomicComposition
getAPOLLOAGORA2Overlap

%% plot various reconstruction statistics for Figure 2 and S1-3
analyseReconstructionStatistics
getMetaboliteStatistics
makeReconstructionSizeBoxplots
plotCorePanReactome
getUniqueReactionsMetabolites

%% Results from the machine learning classifier on strain level
% Extracts the results form the machine learning classifier on the
% strain-level reconstruction data. Creates Table 1 and Table S6.
create_Table_1
create_Table_S6

%% Extraction of properties for microbiome models
getMicrobiomeModelStatistics

%% Results from the machine learning classifier on microbiome level
% Extracts the results from the machine learning classifier on the
% personalised microbiome modeling reconstruction data. Creates Table 2 and
% Table S10. Also exports the top stratifying features to be plotted as a
% heatmap.
create_Table_2
create_Table_S10
exportFeatureDataForHeatmap

% Afterwards, the top stratifying features can be plotted using the R
% script "plotTopFeaturesMicrobiome.R".

%% Plotting of statistical analysis results for the 11 microbiome scenarios
% Extracts the results from statistical analyses and creates a summary
% table and panels for Figure 6, and input for figures.
createTable_StatisticalAnalysisResults
extractSignificantData
plotSignificantSubsystems
create_Table_S11

% Afterwards, the most significant reaction abundances can be plotted using the R
% script "plotSignificantReactionsMicrobiome.R".

%% Plotting of statistically significant microbiome fluxes
% Plots statistically significant metabolite fluxes predicted for
% microbiome scenarios for Figure 7.
plotMicrobiomeFluxes
