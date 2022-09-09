
%% Run the creation and refinement of AGORA2 and all subsequent analyses step by step

% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
% solverOK=changeCobraSolver('gurobi','LP');

rootDir = pwd;
addpath(genpath([rootDir filesep 'Scripts']))
addpath(genpath([rootDir filesep 'input']))

%% Define input variables

infoFilePath = 'AGORA2_infoFile.xlsx';
spreadsheetFolder = [rootDir filesep 'spreadsheets'];
numWorkers = 12;
reconVersion = 'AGORA2';
draftFolder = [rootDir filesep 'draftReconstructions'];
inputDataFolder = [rootDir filesep 'InputData'];
refinedFolder = [rootDir filesep 'refinedReconstructions'];
translatedDraftsFolder = [rootDir filesep 'translatedDraftReconstructions'];
testResultsFolder = [rootDir filesep 'TestResults'];

%% Creaction, refinement, testing, and debugging of AGORA2 reconstructions

% prepare input data
prepareInputData(infoFilePath,'spreadsheetFolder',spreadsheetFolder);

% start the refinement pipeline
% Please ensure the Java Heap Memory is at least 4MB. You may have Out of
% Memory errors otherwise. You can restart the DEMETER pipeline in any case,
% previous progress is being saved.
runPipeline(draftFolder, 'infoFilePath', infoFilePath, 'inputDataFolder', inputDataFolder, 'reconVersion', reconVersion, 'numWorkers', numWorkers);

% translate original AGORA draft reconstructions to enable testing
createMatfileKBaseDraftModels

% run test suite
runTestSuiteTools(refinedFolder, infoFilePath, inputDataFolder, reconVersion, 'numWorkers', numWorkers,'translatedDraftsFolder',translatedDraftsFolder);

% run debugging suite
[debuggingReport, fixedModels, failedModels]=runDebuggingTools(refinedFolder,testResultsFolder,inputDataFolder,infoFilePath,reconVersion,'numWorkers',numWorkers);

% create SBML files for AGORA2 reconstructions, note: very time-consuming
createSBMLFiles

%%
% The following steps can be run from the downloaded AGORA2
% reconstructions, however, the input variables (see above) still need to
% be set.

%% Compute and plot properties of AGORA2 models

% get overview of AGORA2 statistics for Figure 1
getModelStatistics
summarizeReconFeatures
getTaxonStats

% get summary of curation efforts
calculateGenesReactionsRemovedAdded

% compute and plot reaction presences for all AGORA2 reconstructions for Figure 2
propertiesFolder = [rootDir filesep 'modelProperties'];
mkdir(propertiesFolder)

% refined reconstructions
getReactionMetabolitePresence(refinedFolder,propertiesFolder,'AGORA2_refined',numWorkers)
producetSNEPlots(propertiesFolder,infoFilePath,'AGORA2_refined')
computeUptakeSecretion(refinedFolder,propertiesFolder,'AGORA2_refined',{},numWorkers)

% compute and plot reaction presences for subsets of AGORA2 reconstructions
% for Figure 2
clusterClassSubsets

%% Perform validation against three independent experimental datasets
% map other reconstruction resources to AGORA2 nomenclature
mapResourcesToAGORA2Namespace

% then run the comparison
runComparisonOtherResources
computeATPFluxConsistency
plotATPFluxConsistency
exportResults
generateSummaryTable

%% Analyze strain-level drug metabolism capabilities
% analyze comparative genomics data
getGenePresence
getDrugGenePresence
getAddedDrugReactionAverages

% compute and plot single-strain drug fluxes
validateDrugFluxesAgainstExperimentalData
computeDrugConversion
computeDrugYields
plotDrugYields

%% Analyze microbiome-level drug metabolism capabilities
% perform microbiome model building and simulations
mapTaxa
runMgPipe_CRC
runMicrobiomeAnalysis

% export plot microbiome drug metabolism simulation data
exportDrugFluxesForPlot
plotMicrobiomeDrugFluxes
plotFluxesAgainstReactionAbundances
plotFluxesAgainstShadowPrices
correlateFluxesWithTaxa
