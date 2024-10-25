
% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

% path to microbe models
modPath=[pwd filesep 'Models'];

numWorkers=20;

dietFilePath='AverageEuropeanDiet_modified';

computeProfiles=false;

%% run the body site microbiomes one by one

abunFilePaths={
    [pwd filesep 'abundances_from_metagenomes' filesep 'normalized_nasal_cavity_abundances.csv']  'Nasal_cavity';
    [pwd filesep 'abundances_from_metagenomes' filesep 'normalized_vagina_abundances.csv']  'Vagina';
    [pwd filesep 'abundances_from_metagenomes' filesep 'normalized_skin_abundances.csv']  'Skin';
    };

for i=1:size(abunFilePaths,1)
    abunFilePath=abunFilePaths{i,1};
    
    % path where to save results
    resPath=[pwd filesep abunFilePaths{i,2}];
    
    [init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] = initMgPipe(modPath, abunFilePath, computeProfiles, 'dietFilePath', dietFilePath, 'resPath', resPath, 'numWorkers', numWorkers);
end
