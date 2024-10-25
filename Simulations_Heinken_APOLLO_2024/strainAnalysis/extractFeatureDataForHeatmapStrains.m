%% extracts data on the top features stratifying taxa and plots them as heatmaps
% read computed strain-level data and summarize it on different taxon
% levels

versions = {
%     '150k','150k_properties',['150k' filesep '150k_rf_analysis_new'],'Pasolli_genomes_taxonomy_info.xlsx'
%     '90k','90k_properties',['90k' filesep '90k_rf_analysis_new'],'Almeida_genomes_taxonomy_info.xlsx'
    '150k_90k_combined','150k_90k_combined_text_files','150k_90k_combined','Combined_taxonomy_info.xlsx'
    };

datatypes = {
    'ComputedFluxes/InternalProduction_','internal_production_results_'
    'ReactionMetabolitePresence/ReactionPresence_','reaction_presence_results_'
    'ComputedFluxes/UptakeSecretion_','uptake_secretion_results_'
    };

% taxa = {'Phylum','Class','Order','Family','Genus','Species'};
taxa = {'Phylum'};

mkdir(['Model_properties_analysis' filesep 'Summary_for_figures' filesep 'Feature_heatmaps'])

for i=1:size(versions,1)
    infoFile = readInputTableForPipeline(['input' filesep versions{i,4}]);
    for j=1:size(datatypes,1)
        % read the file with taxon information
        for k=1:length(taxa)
            % read the the main classifying features
            topFeatures = readInputTableForPipeline(['Model_properties_analysis' filesep versions{i,3} filesep datatypes{j,2} versions{i,1} filesep taxa{k} filesep 'feature_importance' filesep 'final_feature_importance' '.csv']);

            % read the data
            fid = fopen(['Model_properties_analysis' filesep versions{i,2} filesep datatypes{j,1} versions{i,1} '_refined.txt']);
            tline = fgetl(fid);
            tlines = cell(0,1);
            while ischar(tline)
                tlines{end+1,1} = tline;
                tline = fgetl(fid);
            end
            fclose(fid);
            % get the data to keep
            headers=strsplit(tlines{1,1},'	');
%             [C,IA]=intersect(headers,topFeatures(2:end,1));
            [C,IA]=intersect(headers,topFeatures(2:31,1));

            % put the data together
            dataTmp = {};
            for l=1:length(tlines)
                l
                tline=strsplit(tlines{l,1},'	');
                try
                dataTmp(end+1,1:length(IA)+1) = tline(vertcat(1,IA));
                end
            end

            if j==1
                for l=2:size(dataTmp,1)
                    for m=2:size(dataTmp,2)
                        if str2double(dataTmp{l,m}) > 0.0000001
                            dataTmp{l,m} = '1';
                        else
                            dataTmp{l,m} = '0';
                        end
                    end
                end
            elseif j==3
                dataTmp(1,:) = strrep(dataTmp(1,:),'EX_','');
                dataTmp(1,:) = strrep(dataTmp(1,:),'(e)','');
                for l=2:size(dataTmp,1)
                    for m=2:size(dataTmp,2)
                        if contains(dataTmp{1,m},'secretion')
                            if str2double(dataTmp{l,m}) > 0.0000001
                                dataTmp{l,m} = '1';
                            else
                                dataTmp{l,m} = '0';
                            end
                        elseif contains(dataTmp{1,m},'uptake')
                            if str2double(dataTmp{l,m}) <- 0.0000001
                                dataTmp{l,m} = '-1';
                            else
                                dataTmp{l,m} = '0';
                            end
                        end
                    end
                end
            end
% 
%             for l=2:size(dataTmp,1)
%                 for m=2:size(dataTmp,2)
%                     dataTmp{l,m} = num2str(dataTmp{l,m});
%                 end
%             end

            % summarize the data by taxon to reduce size
            dataSummarized = dataTmp(1,:);
            taxCol = find(strcmp(infoFile(1,:),taxa{k}));
            getTax = unique(infoFile(2:end,taxCol));
            getTax(strncmp(getTax, 'unclassified',  12)) = [];
            getTax(strcmp(getTax, '')) = [];
            getTax(strcmp(getTax, 'N/A')) = [];
            dataSummarized(2:length(getTax)+1,1) = getTax;

            for l=1:length(getTax)
                l
                findTax = infoFile(find(strcmp(infoFile(:,taxCol),getTax{l})),1);
                [~,IA] = intersect(dataTmp(:,1),findTax);
                for m=2:size(dataTmp,2)
                    vals = str2double(dataTmp(IA,m));
                    dataSummarized{l+1,m} = sum(vals)/length(IA);
                end
            end
            cell2csv(['Model_properties_analysis' filesep 'Summary_for_figures' filesep 'Feature_heatmaps' filesep versions{i,1} '_' datatypes{j,2} taxa{k} '.csv'],dataSummarized)
        end
    end
end
