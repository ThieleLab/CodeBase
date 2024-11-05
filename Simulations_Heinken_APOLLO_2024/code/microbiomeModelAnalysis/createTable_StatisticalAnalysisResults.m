
% extract results for statistical analysis of microbiome model comparisons
% to add to Table 2

clear all
rootDir = pwd;

defineScenarios

stats = {
    '_stat_micr.csv'
    '_stat_rxn.csv'
    '_stat_rxnpr.csv'
    '_stat_subsys.csv'
    };

results = {
    'Organism_abundance_'
    'Reaction_abundance_'
    'Reaction_presence_'
    'Subsystem_abundance_'
    };

table = {'','Strain abundance','Reaction abundance','Reaction presence','Subsystem abundance','Fluxes'};

% extract the number of significant figures
for i=1:length(scenarios)
    for j=1:length(stats)
        statResults = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{i} filesep scenarios{i} stats{j}]);
        statResults(1,:) = [];
        findSig = find(cell2mat(statResults(:,end))<0.05);
        table{i+1,1} = scenarios{i};
        data = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{i} filesep results{j} scenarios{i} '.csv']);
        data(1,:) = [];
        % find nonzero data
        summed = sum(cell2mat(data(:,2:end)),2);
        
        table{i+1,j+1} = [num2str(length(findSig)) '/' num2str(length(find(summed>0.0000001)))];
    end
    if any(strcmp(scenarios{i},{'IBD_vs_healthy','Infants_undernourished_vs_healthy','PD_vs_healthy'}))
        statResults = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{i} filesep scenarios{i} '_stat_flux.csv']);
        statResults(1,:) = [];
        findSig = find(cell2mat(statResults(:,end))<0.05);
        table{i+1,5} = [num2str(length(findSig)) '/' num2str(size(statResults,1))];
    end
end
