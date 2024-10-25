% extract subsets of 150k/90k

% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
% solverOK=changeCobraSolver('gurobi','LP');

% 150k
modPath = [pwd filesep 'refinedReconstructions'];

inputDataFolder=[pwd filesep 'inputFiles' filesep];
addpath(inputDataFolder)

infoFilePath=[inputDataFolder 'SequencedGenomesInfo.txt'];
translatedDraftsFolder=[pwd filesep 'translatedDraftReconstructions' filesep];

% name of the reconstruction project
reconVersion='150k';

% define number of workers for parallel computing
numWorkers=20;

% loop through multiple subsets

subHeader = 'Body Site';

% define the feature for which the subset of reconstructions should be 
% extracted.

subFeatures = {
%     'Stool'
%     'Oral cavity'
%     'Skin'
%     'Airways'
    'Vagina'
%     'Nasal cavity'
    };

for i=1:length(subFeatures)
% define the folder where the results from the extracted subset should be
% saved (optional, default folder will be used otherwise).
subsetFolder = [pwd filesep '150k_' strrep(subFeatures{i},' ','')];

% run the extraction of the subset and its analysis
[extractedSubset,subsetFolder] = extractReconstructionResourceSubset(modPath, infoFilePath, subHeader, subFeatures{i}, subsetFolder);

% test the subset of reconstructions if they are already analysis ready

testResultsFolder=[pwd filesep 'Tests_' strrep(subFeatures{i},' ','') filesep];

% quick biomass and ATP check
notGrowing = plotBiomassTestResults(subsetFolder,'testResultsFolder',testResultsFolder, 'numWorkers', numWorkers, 'reconVersion', reconVersion);
tooHighATP = plotATPTestResults(subsetFolder,'testResultsFolder',testResultsFolder, 'numWorkers', numWorkers, 'reconVersion', reconVersion);

runTestSuiteTools(subsetFolder, 'inputDataFolder', inputDataFolder, 'numWorkers', numWorkers, 'testResultsFolder', testResultsFolder, 'infoFilePath', infoFilePath, 'reconVersion', reconVersion)

% get model properties
computeModelProperties(subsetFolder, 'numWorkers', numWorkers, 'propertiesFolder', propertiesFolder, 'infoFilePath', infoFilePath, 'reconVersion', reconVersion, 'customFeatures', customFeatures)

end
