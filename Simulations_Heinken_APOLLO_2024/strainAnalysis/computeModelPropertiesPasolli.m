% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

% define number of workers for parallel computing
numWorkers=20;

% set the path to refined reconstructions
refinedFolder=[rootDir filesep 'PasolliReconstructions' filesep 'refinedReconstructions'];

% name of the reconstruction project
reconVersion='Pasolli';

cd([rootDir filesep 'PasolliReconstructions'])

% get model properties
computeModelProperties(refinedFolder, infoFilePath, reconVersion, 'numWorkers', numWorkers)

cd(rootDir)
