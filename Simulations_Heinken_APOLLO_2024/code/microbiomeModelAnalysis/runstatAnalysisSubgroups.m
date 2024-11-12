%% Code to run statistical analyses for the all the microbiome scenarios in the APOLLO database

% The path to the APOLLO raw data supplementary materials
dataDir = 'D:\APOLLOIndvModels\data_APOLLO\data\analysis_MicrobiomeModels\Scenarios';

% Names of the directories for each scenario
scenarioFolderNames = {
    'T2D_vs_healthy'
    'PD_vs_healthy'
    'Obesity_vs_normalweight'
    'Infection_vs_healthy'
    'Infection_antibiotics_vs_no_antibiotics'
    'Infants_undernourished_vs_healthy'
    'Infants_premature_vs_healthy'
    'IBD_vs_Healthy'
    'Adults_vs_infants_healthy'
    'Adults_healthy_by_country'
    'Adults_body_sites_healthy'};

% Column header to be used in each scenario's metadata file to 
stratColumnName = {
    'DiseaseName'
    'DiseaseName'
    'DiseaseName'
    'DiseaseName'
    'Antibiotics'
    'DiseaseName'
    'DiseaseName'
    'DiseaseName'
    'AgeGroup'
    'Country'
    'BodySite'};

for i = 4:size(scenarioFolderNames,1)
    % Generate the sceneario specfic folder name
    scenarioDirectory = [dataDir, filesep, scenarioFolderNames{i}];
    
    % Enter the specific folder
    cd(scenarioDirectory);
    
    % Run the statistical analysis
    statAnalysisSubgroups(scenarioFolderNames{i}, stratColumnName{i});
end