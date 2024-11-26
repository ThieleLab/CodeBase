
% create a metadata file and folder with the subset of samples for each
% scenario

metadata=readInputTableForPipeline([rootDir filesep 'input' filesep 'Metadata.csv']);
% get the list of microbiome samples in each scenario
scenarioDefinition=readInputTableForPipeline([rootDir filesep 'input' filesep 'ScenarioDefinition.xlsx']);

defineScenarios

mkdir([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios'])
cd([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios'])

for i=1:length(scenarios)
    mkdir(scenarios{i})
    cd(scenarios{i})
    findCol = find(strcmp(scenarioDefinition(1,:),scenarios{i}));
    members = unique(scenarioDefinition(2:end,findCol));
    metadataRed = metadata;
    [C,I] = setdiff(metadataRed(:,1),members,'stable');
    metadataRed(I(2:end),:)=[];
    cell2csv([scenarios{i} '_samples.csv'],metadataRed);
    cd ..
end

cd(rootDir)
