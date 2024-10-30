
% plot significant reactions in microbiome scenarios by subsystem
% by importing the resulting tables into Circos.

defineScenarios

stats = {
    '_stat_rxn.csv'
    '_stat_rxnpr.csv'
    };

% extract the significant features for each scenario
SigFeatures = {};
for i=1:length(stats)
    for j=1:length(scenarios)
        statResults = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{j} filesep scenarios{j} stats{i}]);
        statResults(1,:) = [];
        findNS = find(cell2mat(statResults(:,end))>=0.05);
        statResults(findNS,:) = [];
        findNS = find(isnan(cell2mat(statResults(:,end))));
        statResults(findNS,:) = [];
        SigFeatures{j,i}(1:size(statResults,1),1) = statResults(:,3);
    end
    feats = {};
    for j=1:length(scenarios)
        feats = union(feats,SigFeatures{j,i});
    end
    table = zeros(length(feats),length(scenarios));
    for j=1:length(scenarios)
        for k=1:length(feats)
            table(k,j)=length(find(strcmp(SigFeatures{j,i},feats{k})));
        end
    end
    % exclude exchange/demand and transport
    [C,I] = intersect(feats,{'Exchange/demand','Transport'});
    feats(I,:) = [];
    table(I,:) = [];

    % export the data
    scenarios(find(sum(table,1)==0)) = [];
    table(:,find(sum(table,1)==0)) = [];
    exportTable = cell(length(feats)+1,length(scenarios)+1);
    exportTable(1,2:end) = scenarios;
    exportTable(2:end,1) = feats;
    exportTable(2:end,2:end) = num2cell(table);
    exportTable{1,1} = 'Feature';
    exportTable(2:end,1) = strrep(exportTable(2:end,1),' ','_');
    exportTable(2:end,1) = strrep(exportTable(2:end,1),',','');
    exportTable(2:end,1) = strrep(exportTable(2:end,1),'/','_');
    
    if i==1
        writetable(cell2table(exportTable),[rootDir filesep 'results' filesep 'microbiomes' filesep 'Summary_for_figures' filesep 'SigFeats_ReactionAbundance'],'writeVariableNames',false,'Delimiter','tab')
    else
        writetable(cell2table(exportTable),[rootDir filesep 'results' filesep 'microbiomes' filesep 'Summary_for_figures' filesep 'SigFeats_ReactionPresence'],'writeVariableNames',false,'Delimiter','tab')
    end
end

