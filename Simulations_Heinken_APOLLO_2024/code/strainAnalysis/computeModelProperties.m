
%% compute model properties for APOLLO

% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

% define number of workers for parallel computing
numWorkers=20;

mkdir([rootDir filesep 'data' filesep 'analysis_ModelProperties'])

%% Pasolli reconstructions
% set the path to refined reconstructions
refinedFolder=[rootDir filesep 'data' filesep 'PasolliReconstructions' filesep 'refinedReconstructions'];

% compute model properties
propertiesFolder = [rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Pasolli'];

getReactionMetabolitePresence(refinedFolder,propertiesFolder,[reconVersion '_refined'],numWorkers)
computeUptakeSecretion(refinedFolder,propertiesFolder,[reconVersion '_refined'],{},numWorkers)
computeInternalMetaboliteProduction(refinedFolder,propertiesFolder,[reconVersion '_refined'],{},numWorkers)

%% Almeida reconstructions
% set the path to refined reconstructions
refinedFolder=[rootDir filesep 'data' filesep 'AlmeidaReconstructions' filesep 'refinedReconstructions'];

% compute model properties
propertiesFolder = [rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Almeida'];

getReactionMetabolitePresence(refinedFolder,propertiesFolder,[reconVersion '_refined'],numWorkers)
computeUptakeSecretion(refinedFolder,propertiesFolder,[reconVersion '_refined'],{},numWorkers)
computeInternalMetaboliteProduction(refinedFolder,propertiesFolder,[reconVersion '_refined'],{},numWorkers)
