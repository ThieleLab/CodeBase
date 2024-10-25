
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
    '_stat_micr.csv'
    '_stat_rxn.csv'
    '_stat_rxnpr.csv'
    '_stat_subsys.csv'
    };

results = {
    'Organism_abundance_'
    'Reaction_abundance_'
    'Reactions_presence_'
    'Subsystem_abundance_'
    };

table = {'','Strain abundance','Reaction abundance','Reaction presence','Subsystem abundance','Fluxes'};

% extract the number of signficant figures
for i=1:length(datasets)
    for j=1:length(stats)
        statResults = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{i} filesep datasets{i} stats{j}]);
        statResults(1,:) = [];
        findSig = find(cell2mat(statResults(:,end))<0.05);
        table{i+1,1} = datasets{i};
        data = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{i} filesep results{j} datasets{i} '.csv']);
        data(1,:) = [];
        % find nonzero data
        summed = sum(cell2mat(data(:,2:end)),2);
        
        table{i+1,j+1} = [num2str(length(findSig)) '/' num2str(length(find(summed>0.0000001)))];
    end
    if i==4 || i==6 || i==10
        statResults = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{i} filesep datasets{i} '_stat_flux.csv']);
        statResults(1,:) = [];
        findSig = find(cell2mat(statResults(:,end))<0.05);
        table{i+1,5} = [num2str(length(findSig)) '/' num2str(size(statResults,1))];
    end
end
