
% define datasets with matched c
datasets={
    'Adults_body_sites_healthy' % nasal cavity, skin, and vagina samples
    % the following are all gut samples
    'Adults_by_country' % all healthy adults with country information
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

combinedDatasets={
    [rootDir filesep 'MicrobiomeModeling' filesep 'AllData' filesep 'ModelStatisticsCombined.csv'],'Model sizes','Number'
    [rootDir filesep 'MicrobiomeModeling' filesep 'AllData' filesep 'normalizedAbundanceCombined.csv'],'Organism abundance','Relative abundance'
    [rootDir filesep 'MicrobiomeModeling' filesep 'AllData' filesep 'ReactionAbundanceCombined.csv'],'Reaction abundance','Relative abundance'
    [rootDir filesep 'MicrobiomeModeling' filesep 'AlldData' filesep 'SubsystemAbundanceCombined.csv'],'Subsystem abundance','Relative abundance'
    [rootDir filesep 'MicrobiomeModeling' filesep 'AllData' filesep 'ReactionPresenceCombined.csv'],'Reactions presence','Presence'
    };

cd([rootDir filesep 'Analysis_microbiome_models' filesep 'Datasets'])

for i=1:length(datasets)
    cd(datasets{i})

    metadata=readInputTableForPipeline([datasets{i} '_samples.csv']);

    for j=1:size(combinedDatasets,1)
        j
        dataTable = readInputTableForPipeline(combinedDatasets{j,1});
        if j==2
            data=dataTable;
        else
            data=dataTable';
        end

        % remove samples not in subset
        [C,IA]=setdiff(data(1,:),metadata(2:end,1),'stable');
        data(:,IA(2:end))=[];

        if j==2
            % remove taxa that are not present in the dataset
            cnt=1;
            delArray=[];
            for k=2:size(data,1)
                if sum(str2double(data(k,2:end))) < 0.000001
                    delArray(cnt)=k;
                    cnt=cnt+1;
                end
            end
            data(delArray,:) = [];
        end

        % save the data for the subset of samples
        cell2csv([strrep(combinedDatasets{j,2}, ' ','_') '_' datasets{i} '.csv'],data)
    end
    cd ..
end
