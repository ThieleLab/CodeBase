
% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

modPath=[pwd filesep 'MicrobiomeModels_AED'];

numWorkers = 40;

% some potentially interesting metabolites
metList={'ac','ppa','but','for','isobut','isoval','lac_D','lac_L','etoh','indole','tma','phe_L','tyr_L','trp_L','leu_L','ile_L','val_L'};

mkdir('Strain_Contributions')
resPath = [pwd filesep 'Strain_Contributions'];

predictMicrobeContributions(modPath, 'resPath', resPath, 'metList', metList, 'numWorkers', numWorkers)

infoFilePath = [pwd filesep 'gut_metagenomic_samples_metadata.csv'];
sampleGroupHeaders = {'Disease name','Age group','Location','Country'};
analyzeMgPipeResults(infoFilePath, resPath, 'sampleGroupHeaders', sampleGroupHeaders);
