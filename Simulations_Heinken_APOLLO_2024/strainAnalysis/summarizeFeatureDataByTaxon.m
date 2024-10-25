% summarize data by taxonomy

%% extracts data on the top features stratifying taxa and plots them as heatmaps
% read computed strain-level data and summarize it on different taxon
% levels

versions = {
    %     '150k','150k_properties','SequencedGenomesInfo.txt',['150k' filesep '150k_rf_analysis_new']
    %     '90k','90k_properties','mags-gut_qs50_checkm_adaptedToPipeline.xlsx',['90k' filesep '90k_rf_analysis_new']
    '150k_90k_combined','150k_90k_combined_text_files','Combined_taxonomy_info.csv',['150k' filesep '150k_rf_analysis_new']
    };

taxa = {'Phylum','Class','Order','Family','Genus','Species'};

for i=1:size(versions,1)
    % summarize data on taxon levels

    files = {
        'ReactionPresence_',['ReactionMetabolitePresence' filesep 'reactionPresence_' versions{i,1} '_refined.txt']
        'UptakeSecretion_',['ComputedFluxes' filesep 'UptakeSecretion_' versions{i,1} '_refined_qualitative.txt']
        'InternalProduction_',['ComputedFluxes' filesep 'InternalProduction_' versions{i,1} '_refined_qualitative.txt']
        };

    % if taxonomic information is available, calculate on different taxon levels

    for f=1:size(files,1)
        data = readInputTableForPipeline([versions{i,1} filesep files{f,2}]);
        taxa={'Species','Genus','Family','Order','Class','Phylum'};

        for t=1:length(taxa)
            dataByTaxon=data(1,:);

            % read the the main classifying features and remove all other
            % features
            topFeatures = readInputTableForPipeline(['Model_properties_analysis' filesep versions{i,4} filesep datatypes{j,2} versions{i,1} filesep taxa{t} filesep 'feature_importance' filesep 'important_feature_importance' '.csv']);
            [C,IA] = setdiff(dataByTaxon(:,1),topFeatures(2:end,1),'stable');
            dataByTaxon(IA(2:end),:) = [];

            % get taxonomical information
            infoFile = readInputTableForPipeline(versions{i,3});

            taxCol=find(strcmp(infoFile(1,:),taxa{t}));
            getTax=unique(infoFile(2:end,taxCol));
            getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
            dataByTaxon(2:size(getTax)+1,1)=getTax;
            for j=1:length(getTax)
                getStrains=infoFile(find(strcmp(infoFile(:,taxCol),getTax{j})),1);
                taxonSummary=[];
                for k=1:length(getStrains)
                    findModel=find(strcmp(data(:,1),getStrains{k}));
                    taxonSummary(k,1:size(data,2)-1)=str2double(data(findModel,2:end));
                end
                for l=2:size(dataByTaxon,2)
                    dataByTaxon{j+1,l}=sum(taxonSummary(:,l-1));
                end
            end
            % delete empty columns
            cSums = [NaN,max(cell2mat(dataByTaxon(2:end,2:end)))];
            dataByTaxon(:,find(cSums<0.00000001))=[];
            % normalize the data to the highest abundance value for each
            % taxon
            for j=2:size(dataByTaxon,1)
                maxAll=max(cell2mat(dataByTaxon(j,2:end)));
                for k=2:size(dataByTaxon,2)
                    dataByTaxon{j,k}=dataByTaxon{j,k}/maxAll;
                end
            end
            for j=2:size(dataByTaxon,1)
                for k=2:size(dataByTaxon,2)
                    if isnan(dataByTaxon{j,k})
                        dataByTaxon{j,k}=0;
                    end
                end
            end
            cell2csv(['Model_properties_analysis' filesep 'Summary_for_figures' filesep 'Feature_heatmaps' filesep versions{i,1} '_' datatypes{j,1} '_' taxa{t} '.csv'],dataByTaxon)
        end
    end
end
