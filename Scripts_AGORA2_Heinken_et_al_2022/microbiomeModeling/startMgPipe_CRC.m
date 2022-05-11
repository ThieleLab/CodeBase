
mkdir('Modeling_CRC')

%% normalize coverage
abunFilePath = 'coverage_CRC.csv';
[normalizedCoverage,normalizedCoveragePath] = normalizeCoverage(abunFilePath,0.01);

%% Create models with Japanese diet
% path to microbe models
modPath=[rootDir filesep 'panModelsAGORA2' filesep 'Species'];

% path where to save results
resPath=[rootDir filesep 'Modeling_CRC' filesep 'MicrobiomeModels'];

% path to and name of the file with dietary information.
dietFilePath='JapaneseDiet';

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, normalizedCoveragePath, false, 'dietFilePath', dietFilePath, 'resPath', resPath, 'numWorkers', numWorkers);

% move diet models to a new path
dInfo = dir([resPath filesep 'Diet']);
modelList={dInfo.name};
modelList=modelList';
modelList(~(contains(modelList(:,1),'.mat')),:)=[];
for i=1:length(modelList)
    movefile([resPath filesep 'Diet' filesep modelList{i}],[rootDir filesep 'Modeling_CRC' filesep 'JapaneseDiet' filesep modelList{i}])
end

%% compute uptake and secretion profiles for Japanese diet

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, normalizedCoveragePath, true, 'dietFilePath', dietFilePath, 'resPath', resPath,'numWorkers', 4);

%% create models with European diet

% path to and name of the file with dietary information.
dietFilePath='AverageEuropeanDiet';

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(modPath, normalizedCoveragePath, false, 'dietFilePath', dietFilePath, 'resPath', resPath, 'numWorkers', numWorkers);

% move diet models to a new path
dInfo = dir([resPath filesep 'Diet']);
modelList={dInfo.name};
modelList=modelList';
modelList(~(contains(modelList(:,1),'.mat')),:)=[];
for i=1:length(modelList)
    movefile([resPath filesep 'Diet' filesep modelList{i}],[rootDir filesep 'Modeling_CRC' filesep 'AverageEuropeanDiet' filesep modelList{i}])
end
