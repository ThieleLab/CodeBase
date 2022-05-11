
% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
% prevent creation of log files
changeCobraSolverParams('LP', 'logFile', 0);

numWorkers=12;

% define folders where data will be retrieved/saved
rootDir = strrep(matlab.desktop.editor.getActiveFilename,'runAGORA2.m','');
draftFolder = [rootDir filesep 'Input_Models'];
refinedFolder = [rootDir filesep 'Output_Models'];
translatedDraftsFolder = [rootDir filesep 'Draft_Reconstructions_Matfiles'];
summaryFolder = [rootDir filesep 'curationSummary'];
testResultsFolder = [rootDir filesep 'TestResults'];
propertiesFolder = [rootDir filesep 'ModelProperties'];

spreadsheetFolder = [rootDir filesep 'spreadsheets'];

inputDataFolder=[rootDir filesep 'InputFiles'];
infoFilePath = 'AGORA2_infoFile.xlsx';

reconVersion='AGORA2';

addpath(genpath([rootDir filesep 'Scripts']))

% prepare input data including comparative genomics data
prepareInputData(infoFilePath,'inputDataFolder',inputDataFolder,'spreadsheetFolder',spreadsheetFolder);

% run refinement of draft reconstructions resulting in AGORA2
runPipeline(draftFolder, 'infoFilePath', infoFilePath, 'inputDataFolder', inputDataFolder, 'refinedFolder', refinedFolder, 'translatedDraftsFolder', translatedDraftsFolder, 'summaryFolder', summaryFolder, 'numWorkers', numWorkers, 'reconVersion', reconVersion)

% run test suite
% convert original AGORA draft reconstructions to namespace
createMatfileKBaseDraftModels
runTestSuiteTools(refinedFolder, infoFilePath, inputDataFolder, reconVersion, 'numWorkers', numWorkers, 'testResultsFolder', testResultsFolder,'translatedDraftsFolder',translatedDraftsFolder)
[debuggingReport, fixedModels, failedModels]=runDebuggingTools(refinedFolder,testResultsFolder,inputDataFolder,infoFilePath,reconVersion,'numWorkers',numWorkers);

% get model properties
computeModelProperties(refinedFolder, infoFilePath, reconVersion, 'translatedDraftsFolder', translatedDraftsFolder, 'numWorkers', numWorkers, 'propertiesFolder', propertiesFolder)

% get reconstruction features-for Figure 1
clusterClassSubsets
getReconstructionStatistics
prepareSummaryTableReconFeatures
getTaxonStats

% get summary of curation efforts
calculateGenesReactionsRemovedAdded
curationStatus = getCurationStatus(infoFilePath,inputDataFolder,0);
getSums = sum(cell2mat(curationStatus(2:end,2:end))==2,2)>=1;
numCurated = sum(getSums);
percentCurated = sum(getSums)/length(getSums);

% compare with other resources against experimental data
% map reaction identifiers to AGORA2 namespace to enable comparison
mapResourcesToAGORA2Namespace

% extract the data
cd([rootDir filesep 'Comparison_other_GEMs' filesep 'inputData'])
convertLimDataToTable
convertMadinDataToTable
runMappingBacDiveData
addpath(genpath([rootDir filesep 'Comparison_other_GEMs' filesep 'inputData']))
cd(rootDir)

% then run the comparison
runComparisonOtherResources
cd(rootDir)
computeATPFluxConsistency
plotATPFluxConsistency
exportResults
generateSummaryTable

% get some AGORA2 reconstructions for MEMOTE
getRandomRecons

% analyze comparative genomics data
getGenePresence
getDrugGenePresence
getAddedDrugReactionAverages

% compute and plot single-strain drug fluxes
validateDrugFluxesAgainstExperimentalData
computeDrugConversion
computeDrugYields
plotDrugYields

% perform microbiome modeling
calculateMicrobiomeCoverage
makePanModels
startMgPipe_CRC
startMicrobiomeAnalysis
correlateFluxesWithTaxa

% plot microbiome simulation data
exportDrugFluxesForPlot
plotDrugViolins
plotFluxesAgainstReactionAbundances
plotFluxesAgainstShadowPrices
