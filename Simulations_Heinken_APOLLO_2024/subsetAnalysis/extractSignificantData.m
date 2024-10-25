
% plot significant reactions in microbiome datasets by subsystem

mkdir(['Analysis_microbiome_models' filesep 'Summary_for_figures' filesep 'SignificantData'])

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
    for j=1:length(datasets)
        statResults = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{j} filesep datasets{j} stats{i}]);
        statResults(1,:) = [];
        findNS = find(cell2mat(statResults(:,end))>=0.05);
        statResults(findNS,:) = [];
        findNS = find(isnan(cell2mat(statResults(:,end))));
        statResults(findNS,:) = [];
        data = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{j} filesep results{i} datasets{j} '.csv']);
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
        writetable(cell2table(data),['Analysis_microbiome_models' filesep 'Summary_for_figures' filesep 'SignificantData' filesep results{i} datasets{j} '.csv'],'writeVariableNames',false)
    end
end

