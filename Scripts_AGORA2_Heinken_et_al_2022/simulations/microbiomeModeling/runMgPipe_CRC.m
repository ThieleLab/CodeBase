
mkdir('Modeling_CRC')

%% create pan-species models

% path to save pan-species models
spPath=[rootDir filesep 'Modeling_CRC' filesep 'Species'];
createPanModels(refinedFolder,spPath,'Species',numWorkers,'AGORA2_infoFile.xlsx');

% test biomass and ATP production in pan-species models
testResultsFolder=[rootDir filesep 'Modeling_CRC' filesep 'TestResults_Species'];
[notGrowing,biomassFluxes] = plotBiomassTestResults(spPath, 'Species', 'numWorkers',numWorkers, 'testResultsFolder', testResultsFolder);
[tooHighATP,atpFluxes] = plotATPTestResults(spPath, 'Species', 'numWorkers',numWorkers, 'testResultsFolder', testResultsFolder);

%% Create models with Japanese diet

% path to file with normalized abundances
abunFilePath = [rootDir filesep 'normalizedCoverage.csv'];

% path where to save results
resPath=[rootDir filesep 'Modeling_CRC' filesep 'MicrobiomeModels'];

% path to and name of the file with dietary information.
dietFilePath='JapaneseDiet';

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(spPath, abunFilePath, false, 'dietFilePath', dietFilePath, 'resPath', resPath, 'numWorkers', numWorkers);

% move diet-constrained models to a new path
mkdir([rootDir filesep 'Modeling_CRC' filesep 'JapaneseDiet'])
dInfo = dir([resPath filesep 'Diet']);
modelList={dInfo.name};
modelList=modelList';
modelList(~(contains(modelList(:,1),'.mat')),:)=[];
for i=1:length(modelList)
    movefile([resPath filesep 'Diet' filesep modelList{i}],[rootDir filesep 'Modeling_CRC' filesep 'JapaneseDiet' filesep modelList{i}])
end

%% create models with European diet

% path to and name of the file with dietary information.
dietFilePath='AverageEuropeanDiet';

[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(spPath, normalizedCoveragePath, false, 'dietFilePath', dietFilePath, 'resPath', resPath, 'numWorkers', numWorkers);

% move diet-constrained models to a new path
mkdir([rootDir filesep 'Modeling_CRC' filesep 'AverageEuropeanDiet'])
dInfo = dir([resPath filesep 'Diet']);
modelList={dInfo.name};
modelList=modelList';
modelList(~(contains(modelList(:,1),'.mat')),:)=[];
for i=1:length(modelList)
    movefile([resPath filesep 'Diet' filesep modelList{i}],[rootDir filesep 'Modeling_CRC' filesep 'AverageEuropeanDiet' filesep modelList{i}])
end

%% compute complete uptake and secretion profiles of the microbiome models on Japanese diet
% Warning: this step is time-consuming.

poolobj = gcp('nocreate');
delete(poolobj)

% path to pan-species models
spPath=[rootDir filesep 'Modeling_CRC' filesep 'Species'];

% path to file with normalized abundances
abunFilePath = [rootDir filesep 'normalizedCoverage.csv'];

% path where to save results
resPath=[rootDir filesep 'Modeling_CRC' filesep 'MicrobiomeModels'];

dietFilePath='JapaneseDiet';
[init, netSecretionFluxes, netUptakeFluxes, Y] = initMgPipe(spPath, abunFilePath, true, 'dietFilePath', dietFilePath, 'resPath', resPath,'numWorkers', 4);

