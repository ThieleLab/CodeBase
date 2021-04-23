% This scripts creates personalized microbiome models constrained with a 
% simulated Average European diet for the COMBO/PLEASE cohort. Note that
% due to changes in the diet implementation in the Microbiome Modeling
% Toolbox, they may slightly differ from the models originally created for
% the paper.

% Almut Heinken, 03/21

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

%% Set input variables

% Do not perform analysis in mgPipe in this case, just export the models
% with diet constraints
computeProfiles = false;

% path where to save constrained models
resPath = [rootDir filesep 'MicrobiomeModels'];

% path to and name of the file with dietary information
dietFilePath = 'AverageEuropeanDiet';

% save models with diet constrains implemented
saveConstrModels = true;

%% Run model building pipeline

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(agoraPath, abunFilePath, computeProfiles, 'resPath', resPath, 'dietFilePath', dietFilePath, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers);
