
% extract the results from statistical analysis for microbiome models for
% Table S11a-d

clear all
rootDir = pwd;

mkdir([rootDir filesep 'results' filesep 'microbiomes' filesep 'Statistical_results'])

defineScenarios

% determine the subset for which flux data exists
fluxSubsets={
    'IBD_vs_healthy' % from PMID:24629344
    'Infants_undernourished_vs_healthy' % undernourished and normal infants from Bangladesh
    'PD_vs_healthy' % from PMID:28662719
    };

stats = {
    '_stat_rxn.csv'
    '_stat_rxnpr.csv'
    '_stat_subsys.csv'
    '_stat_flux.csv'
    };

mkdir([rootDir filesep 'results' filesep 'microbiomes' filesep 'Statistical_results'])

% create a summary table
for i=1:length(stats)
    if i<4
        feats = {};
        for j=1:length(scenarios)
            statResults = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{j} filesep scenarios{j} stats{i}]);
            feats = union(feats,statResults(2:end,1));
        end
        if i<3
            table = {'Reaction','Subsystem','Subsystem general'};
            cnt=3;
        else
            table = {'Subsystem'};
            cnt=1;
        end
        table(2:length(feats)+1,1)=feats;

        for j=1:length(scenarios)
            statResults = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{j} filesep scenarios{j} stats{i}]);
            table{1,cnt+1} = ['Skewness ' scenarios{j}];
            table{1,cnt+2} = ['Kurtosis ' scenarios{j}];
            table{1,cnt+3} = ['p-value ' scenarios{j}];
            table{1,cnt+4} = ['After FDR correction ' scenarios{j}];
            for k=2:size(statResults,1)
                findRxn = find(strcmp(table(:,1),statResults{k,1}));
                table{findRxn,2} = statResults{k,2};
                table{findRxn,3} = statResults{k,3};
                table{findRxn,cnt+3} = statResults{k,size(statResults,2)-1};
                table{findRxn,cnt+4} = statResults{k,size(statResults,2)};
                table{findRxn,cnt+1} = statResults{k,4};
                table{findRxn,cnt+2} = statResults{k,5};
            end
            cnt = cnt + 4;
        end
    else
        feats = {};
        for j=1:length(fluxSubsets)
            statResults = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep fluxSubsets{j} filesep fluxSubsets{j} stats{i}]);
            feats = union(feats,statResults(2:end,1));
        end

        table = {'Flux'};
        cnt=1;
        table(2:length(feats)+1,1)=feats;

        for j=1:length(fluxSubsets)
            statResults = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep fluxSubsets{j} filesep fluxSubsets{j} stats{i}]);
            table{1,cnt+1} = ['Skewness ' fluxSubsets{j}];
            table{1,cnt+2} = ['Kurtosis ' fluxSubsets{j}];
            table{1,cnt+3} = ['p-value ' fluxSubsets{j}];
            table{1,cnt+4} = ['After FDR correction ' fluxSubsets{j}];
            for k=2:size(statResults,1)
                findRxn = find(strcmp(table(:,1),statResults{k,1}));
                table{findRxn,2} = statResults{k,2};
                table{findRxn,3} = statResults{k,3};
                table{findRxn,cnt+3} = statResults{k,size(statResults,2)-1};
                table{findRxn,cnt+4} = statResults{k,size(statResults,2)};
                table{findRxn,cnt+1} = statResults{k,4};
                table{findRxn,cnt+2} = statResults{k,5};
            end
            cnt = cnt + 4;
        end
    end
   writetable(cell2table(table),[rootDir filesep 'results' filesep 'microbiomes' filesep 'Statistical_results' filesep 'Summary_' stats{i} '.csv'],'writeVariableNames',false)
end
