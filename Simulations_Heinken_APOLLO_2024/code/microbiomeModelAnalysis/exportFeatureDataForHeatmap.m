% get the top hits from random forests analysis and plot their distribution
% in the samples

clear all
rootDir = pwd;

mkdir([rootDir filesep 'results' filesep 'microbiomes' filesep 'Summary_for_figures' filesep 'Feature_heatmaps'])

% define the four types of data
datatypes={
    'reaction_abundance','Reaction_abundance_'
    'reaction_presence','Reactions_presence_'
    'subsystem_abundance','Subsystem_abundance_'
    };

defineScenarios

for i=1:length(scenarios)
    for j=1:size(datatypes,1)
        % read the data and the main classifying features
        data = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' scenarios{i} filesep datatypes{j,2} scenarios{i} '.csv']);
        topFeatures = readInputTableForPipeline([rootDir filesep 'data' filesep 'RF_Results_Microbiome' filesep scenarios{i} filesep 'feature_importance' filesep 'final_feature_importance_' datatypes{j,1} '.csv']);
        [C,IA] = setdiff(data(:,1),topFeatures(2:end,1),'stable');
        data(IA(2:end),:) = [];
        cell2csv([rootDir filesep 'results' filesep 'microbiomes' filesep 'Summary_for_figures' filesep 'Feature_heatmaps' filesep datatypes{j,2} scenarios{i} '.csv'],data)
    end
end
