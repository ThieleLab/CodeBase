
% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

% path to microbe models
modPath=[rootDir filesep 'data' filesep 'BodySiteMicrobiomes' filesep 'Models'];

numWorkers=20;

dietFilePath='AverageEuropeanDiet_modified';

computeProfiles=false;

%% run the body site microbiomes one by one

abunFilePaths={
    [rootDir filesep 'input' filesep 'normalized_nasal_cavity_abundances.csv']  'Nasal_cavity';
    [rootDir filesep 'input' filesep 'normalized_vagina_abundances.csv']  'Vagina';
    [rootDir filesep 'input' filesep 'normalized_skin_abundances.csv']  'Skin';
    };

for i=1:size(abunFilePaths,1)
    abunFilePath=abunFilePaths{i,1};
    
    % path where to save results
    resPath=[rootDir filesep 'data' filesep 'BodySiteMicrobiomes' filesep abunFilePaths{i,2}];
    
    [init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] = initMgPipe(modPath, abunFilePath, computeProfiles, 'dietFilePath', dietFilePath, 'resPath', resPath, 'numWorkers', numWorkers);
end
