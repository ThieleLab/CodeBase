
% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

% path to microbe models
modPath='D:\90000_genomes_Almeida_2019\Models';

numWorkers=4;

computeProfiles=true;

resPath = [pwd filesep 'MicrobiomeModels'];

abunFilePath=[pwd filesep 'Organism_abundance_PD_vs_healthy.csv'];

infoFilePath = [pwd filesep 'PD_vs_healthy_samples.csv'];

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] = initMgPipe(modPath, abunFilePath, computeProfiles, 'resPath', resPath, 'numWorkers', numWorkers);

% create models with diet constraints
implementDiet

% compute selected metabolites of interest
computeShadowPrices

