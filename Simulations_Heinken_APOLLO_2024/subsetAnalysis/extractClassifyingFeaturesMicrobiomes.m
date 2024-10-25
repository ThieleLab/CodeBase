
% extract classifying features for Table S10
table = {};
cnt=1;

% define the four types of data
datatypes={
    'Organism abundance_','organism_abundance','Organism_abundance_'
   'Reaction abundance_','reaction_abundance','Reaction_abundance_'
    'Reaction presence','reaction_presence','Reactions_presence_'
    'Subsystem abundance','subsystem_abundance','Subsystem_abundance_'
    };

% define the different dataset results to extract
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

for i=1:length(datasets)
    for j=1:size(datatypes,1)
        data = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{i} filesep datatypes{j,3} datasets{i} '.csv']);
        metadata = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{i} filesep datasets{i} '_samples.csv']);
        if strcmp(datasets{i},'Infection_antibiotics_vs_no_antibiotics')
            strat='Antibiotics';
        elseif strcmp(datasets{i},'Infection_resistant_vs_susceptible')
            strat='Stratification';
        elseif strcmp(datasets{i},'Adults_body_sites_healthy')
            strat='Body site';
        elseif strcmp(datasets{i},'Adults_vs_infants_healthy')
            strat='Age group';
        else
            strat='Disease name';
        end

        % read the data and the main classifying features
        feats = table2cell(readtable(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'RF_Results' filesep datasets{i} filesep 'feature_importance' filesep 'final_feature_importance_' datatypes{j,2} '.csv']));
        if size(feats,2)>2
            % work around for tables not read properly
            for k=1:size(feats,1)
                feats{k,2}=num2str(feats{k,2});
                [a]=strsplit(feats{k,size(feats,2)},',');
                feats{k,size(feats,2)-1}=a{1};
                feats{k,size(feats,2)}=a{2};
                feat='';
                for l=1:size(feats,2)-2
                feat=[feat feats{k,l} '_'];
                end
                feat=[feat feats{k,size(feats,2)-1}];
                feats{k,1}=feat;   
            end
            feats(:,2:size(feats,2)-1)=[];
        end
        table{1,cnt} = [datasets{i,1} ' ' lower(datatypes{j,1})];
        table(2:size(feats,1)+1,cnt) = feats(:,1);
        table(2:size(feats,1)+1,cnt+1) = feats(:,2);
        cnt=cnt+2;
    end
end

writetable(cell2table(table),'Table_S10_Feature_Summary.csv','WriteVariableNames',false);

