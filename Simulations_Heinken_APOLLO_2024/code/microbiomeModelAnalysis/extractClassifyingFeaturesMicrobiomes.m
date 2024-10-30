
% extract classifying features for microbiomes for Table 2
table = {};
cnt=1;

% define the four types of data
datatypes={
    'Organism abundance_','organism_abundance','Organism_abundance_'
    'Reaction abundance_','reaction_abundance','Reaction_abundance_'
    'Reaction presence','reaction_presence','Reactions_presence_'
    'Subsystem abundance','subsystem_abundance','Subsystem_abundance_'
    };

defineScenarios

for i=1:length(scenarios)
    for j=1:size(datatypes,1)
        data = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{i} filesep datatypes{j,3} scenarios{i} '.csv']);
        metadata = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios' filesep scenarios{i} filesep scenarios{i} '_samples.csv']);
        if strcmp(scenarios{i},'Infection_antibiotics_vs_no_antibiotics')
            strat='Antibiotics';
        elseif strcmp(scenarios{i},'Infection_resistant_vs_susceptible')
            strat='Stratification';
        elseif strcmp(scenarios{i},'Adults_body_sites_healthy')
            strat='Body site';
        elseif strcmp(scenarios{i},'Adults_vs_infants_healthy')
            strat='Age group';
        else
            strat='Disease name';
        end
        
        % read the data and the main classifying features
        feats = table2cell(readtable([rootDir filesep 'data' filesep 'RF_Results_Microbiome' filesep scenarios{i} filesep 'feature_importance' filesep 'final_feature_importance_' datatypes{j,2} '.csv']));
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
        table{1,cnt} = [scenarios{i,1} ' ' lower(datatypes{j,1})];
        table(2:size(feats,1)+1,cnt) = feats(:,1);
        table(2:size(feats,1)+1,cnt+1) = feats(:,2);
        cnt=cnt+2;
    end
end

writetable(cell2table(table),[rootDir filesep 'results' filesep 'microbiomes' filesep 'Table_2_Feature_Summary.csv'],'WriteVariableNames',false);
