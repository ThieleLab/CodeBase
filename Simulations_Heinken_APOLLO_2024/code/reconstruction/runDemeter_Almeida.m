
clear all
rootDir = pwd;

% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

% define number of workers for parallel computing
numWorkers=20;

mkdir([rootDir filesep 'data' filesep 'AlmeidaReconstructions'])

% set the path to draft reconstructions
draftFolder=[rootDir filesep 'draftReconstructions' filesep 'Almeida' filesep];

% name of the reconstruction project
reconVersion='Almeida';

% prepare input data
infoFilePath=[rootDir filesep 'input' filesep 'Almeida_genomes_taxonomy_info.txt'];
[adaptedInfoFilePath,inputDataFolder] = prepareInputData(infoFilePath);

refinedFolder = [rootDir filesep 'data' filesep 'AlmeidaReconstructions' filesep 'refinedReconstructions'];

% run DEMETER
runDemeter(draftFolder, 'infoFilePath', infoFilePath, 'numWorkers', numWorkers, 'refinedFolder', refinedFolder, 'reconVersion', reconVersion, 'createSBML', true);

% run test suite
[testResultsFolder,curationReport] = runTestSuiteTools(refinedFolder, infoFilePath, inputDataFolder, reconVersion, 'numWorkers', numWorkers);

% run debugging suite
[debuggingReport, fixedModels, failedModels]=runDebuggingTools(refinedFolder,testResultsFolder,infoFilePath,inputDataFolder,reconVersion,'numWorkers',numWorkers);
