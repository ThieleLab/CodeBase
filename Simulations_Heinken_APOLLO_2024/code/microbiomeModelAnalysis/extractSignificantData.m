
% plot significant reactions in microbiome scenarios by subsystem

mkdir([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Summary_for_figures'])
mkdir([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Summary_for_figures' filesep 'SignificantData'])

defineScenarios

stats = {
    '_stat_rxn.csv'
    '_stat_rxnpr.csv'
    '_stat_subsys.csv'
    };

results = {
    'Reaction_abundance_'
    'Reactions_presence_'
    'Subsystem_abundance_'
    };

% extract the significant features for each dataset
SigFeatures = {};
for i=1:length(stats)
    for j=1:length(scenarios)
        statResults = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{j} filesep scenarios{j} stats{i}]);
        statResults(1,:) = [];
        findNS = find(cell2mat(statResults(:,end))>=0.05);
        statResults(findNS,:) = [];
        findNS = find(isnan(cell2mat(statResults(:,end))));
        statResults(findNS,:) = [];
        data = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{j} filesep results{i} scenarios{j} '.csv']);
        if i==3
            % workaround for changed names
            oldNames=data(:,1);
            data(:,1) = strrep(data(:,1),' ','');
            data(:,1) = upper(data(:,1));
            statResults(:,1) = upper(statResults(:,1));
        end
        [C,I] = setdiff(data(:,1),statResults(:,1),'stable');
        if i==3
            data(:,1) = oldNames;
        end
        data(I(2:end),:) = [];
        delArray = [];
        cnt=1;
        for k=2:size(data,1)
            if sum(cell2mat(data(k,2:end)))<0.0000001
                delArray(cnt,1)=k;
                cnt=cnt+1;
            end
        end
        data(delArray,:)=[];
        writetable(cell2table(data),[rootDir filesep 'results' filesep 'microbiomes' filesep 'Summary_for_figures' filesep 'SignificantData' filesep results{i} scenarios{j} '.csv'],'writeVariableNames',false)
    end
end

