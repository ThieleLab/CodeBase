
% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

% define number of workers for parallel computing
numWorkers=20;

rootDir = pwd;

mkdir([rootDir filesep 'PasolliReconstructions'])
cd([rootDir filesep 'PasolliReconstructions'])

% set the path to draft reconstructions
draftFolder=[rootDir filesep 'draftReconstructions' filesep 'Pasolli' filesep];

% name of the reconstruction project
reconVersion='Pasolli';

% prepare input data
infoFilePath=[rootDir filesep 'input' filesep 'Pasolli_genomes_taxonomy_info.xlsx'];
[adaptedInfoFilePath,inputDataFolder] = prepareInputData(infoFilePath);

% run pipeline
runPipeline(draftFolder, 'infoFilePath', adaptedInfoFilePath, 'numWorkers', numWorkers, 'reconVersion', reconVersion);

% run test suite
[testResultsFolder,curationReport] = runTestSuiteTools(refinedFolder, infoFilePath, inputDataFolder, reconVersion, 'numWorkers', numWorkers);

% run debugging suite
[debuggingReport, fixedModels, failedModels]=runDebuggingTools(refinedFolder,testResultsFolder,infoFilePath,inputDataFolder,reconVersion,'numWorkers',numWorkers);

% get model properties
computeModelProperties(refinedFolder, infoFilePath, reconVersion, 'numWorkers', numWorkers)

cd(rootDir)
