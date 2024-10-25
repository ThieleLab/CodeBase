
% plot significant reactions in microbiome datasets by subsystem
% by importing the resulting tables into Circos.ca

datasets={
    'Adults_body_sites_healthy' % nasal cavity, skin, and vagina samples
    'Adults_healthy_by_country' % healthy adult gut samples by country
    'Adults_vs_infants_healthy' % all healthy adults and infants
    'IBD_vs_healthy' % from PMID:24629344
    'Infants_premature_vs_healthy' % all healthy and premature infants
    'Infants_undernourished_vs_healthy' % undernourished and normal infants from Bangladesh
    'Infection_antibiotics_vs_no_antibiotics' % Cholera study, no REF
    'Infection_vs_healthy' % all healthy gut samples vs. infection
    'Obesity_vs_normalweight' % from PMID:23985870, other samples with BMI available
    'PD_vs_healthy' % from PMID:28662719
    'T2D_vs_healthy' % all T2D vs healthy adults
    };

stats = {
    '_stat_rxn.csv'
    '_stat_rxnpr.csv'
    };

% extract the significant features for each dataset
SigFeatures = {};
for i=1:length(stats)
    for j=1:length(datasets)
        statResults = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{j} filesep datasets{j} stats{i}]);
        statResults(1,:) = [];
        findNS = find(cell2mat(statResults(:,end))>=0.05);
        statResults(findNS,:) = [];
        findNS = find(isnan(cell2mat(statResults(:,end))));
        statResults(findNS,:) = [];
        SigFeatures{j,i}(1:size(statResults,1),1) = statResults(:,3);
    end
    feats = {};
    for j=1:length(datasets)
        feats = union(feats,SigFeatures{j,i});
    end
    table = zeros(length(feats),length(datasets));
    for j=1:length(datasets)
        for k=1:length(feats)
            table(k,j)=length(find(strcmp(SigFeatures{j,i},feats{k})));
        end
    end
    % exclude exchange/demand and transport
    [C,I] = intersect(feats,{'Exchange/demand','Transport'});
    feats(I,:) = [];
    table(I,:) = [];

    % export the data
    datasets(find(sum(table,1)==0)) = [];
    table(:,find(sum(table,1)==0)) = [];
    exportTable = cell(length(feats)+1,length(datasets)+1);
    exportTable(1,2:end) = datasets;
    exportTable(2:end,1) = feats;
    exportTable(2:end,2:end) = num2cell(table);
    exportTable{1,1} = 'Feature';
    exportTable(2:end,1) = strrep(exportTable(2:end,1),' ','_');
    exportTable(2:end,1) = strrep(exportTable(2:end,1),',','');
    exportTable(2:end,1) = strrep(exportTable(2:end,1),'/','_');
    
    if i==1
        writetable(cell2table(exportTable),['Analysis_microbiome_models' filesep 'Summary_for_figures' filesep 'SigFeats_ReactionAbundance'],'writeVariableNames',false,'Delimiter','tab')
    else
        writetable(cell2table(exportTable),['Analysis_microbiome_models' filesep 'Summary_for_figures' filesep 'SigFeats_ReactionPresence'],'writeVariableNames',false,'Delimiter','tab')
    end
end

