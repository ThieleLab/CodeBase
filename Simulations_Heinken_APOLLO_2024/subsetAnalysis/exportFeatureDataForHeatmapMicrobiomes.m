% get the top hits from random forests analysis and plot their distribution
% in the samples

mkdir(['Analysis_microbiome_models' filesep 'Summary_for_figures' filesep 'Feature_heatmaps'])

% define the four types of data
datatypes={
    'reaction_abundance','Reaction_abundance_'
    'reaction_presence','Reactions_presence_'
    'subsystem_abundance','Subsystem_abundance_'
    };

% define the different subgroup results to extract
datasets={
    'Adults_body_sites_healthy' % nasal cavity, skin, and vagina samples
    'Adults_healthy_by_country'
    'Adults_vs_infants_healthy' % all healthy adults and infants
    'IBD_vs_healthy' % from PMID:24629344
    'Infants_premature_vs_healthy' % all healthy and premature infants
    'Infants_undernourished_vs_healthy' % undernourished and normal infants from Bangladesh
    'Infection_antibiotics_vs_no_antibiotics' % Cholera study, no REF
    'Infection_resistant_vs_susceptible' % from PMID:30057943
    'Infection_vs_healthy' % all healthy gut samples vs. infection
    'Obesity_vs_normalweight' % from PMID:23985870, other samples with BMI available
    'PD_vs_healthy' % from PMID:28662719
    'T2D_vs_healthy' % all T2D vs healthy adults
    };

for i=1:length(datasets)
    for j=1:size(datatypes,1)
        % read the data and the main classifying features
        data = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{i} filesep datatypes{j,2} datasets{i} '.csv']);
        topFeatures = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'RF_Results' filesep datasets{i} filesep 'feature_importance' filesep 'final_feature_importance_' datatypes{j,1} '.csv']);
%         [C,IA] = setdiff(data(:,1),topFeatures(2:end,1),'stable');
%         if size(topFeatures,1)>30
%             [C,IA] = setdiff(data(:,1),topFeatures(2:31,1),'stable');
%         else
            [C,IA] = setdiff(data(:,1),topFeatures(2:end,1),'stable');
%         end
        data(IA(2:end),:) = [];
        cell2csv(['Analysis_microbiome_models' filesep 'Summary_for_figures' filesep 'Feature_heatmaps' filesep datatypes{j,2} datasets{i} '.csv'],data)
    end
end
