
% read subsystem abundances and summarize by taxon

subs_90k = readInputTableForPipeline('SubsystemPresence_90k_refined.txt');

subs_150k = readInputTableForPipeline('SubsystemPresence_150k_refined.txt');

% rearrange columns
subs_90k{1,1} = '';
subs_150k{1,1} = '';
subs1 = subs_90k(1,:);
subs2 = subs_150k(1,:);
% remove subsystems only in one dataset
[C,IA] = setdiff(subs1,subs2);
subs_90k(:,IA) = [];
[C,IA] = setdiff(subs2,subs1);
subs_150k(:,IA) = [];
subs_all = vertcat(subs_150k,subs_90k(2:end,:));
writetable(cell2table(subs_all),'SubsystemPresence_150k_all_refined.txt','FileType','text','WriteVariableNames',false,'Delimiter','tab');

mkdir('SubsystemsByTaxon')

% define paths with taxon information
paths = {
%     'Pasolli','SequencedGenomesInfo.txt','SubsystemPresence_150k_refined.txt'
%     'Almeida','mags-gut_qs50_checkm_adaptedToPipeline.xlsx','SubsystemPresence_90k_refined.txt'
    'GlobalBiome','Combined_taxonomy_info.xlsx','SubsystemPresence_all_refined.txt'
    };

taxLevels={'Species','Genus','Family','Order','Class','Phylum'};

for p = 1:size(paths,1)
    infoFile = readInputTableForPipeline(paths{p,2});
    data = readInputTableForPipeline(paths{p,3});

    for i=1:length(taxLevels)
        dataByTaxon=data(1,:);
        taxCol=find(strcmp(infoFile(1,:),taxLevels{i}));
        getTax=unique(infoFile(2:end,taxCol));
        if i==1
            getTax(find(strncmp(getTax,'Candidatus',length('Candidatus'))),:)=[];
        end
        getTax(find(strcmp(getTax,'N/A')),:)=[];
        getTax(find(contains(getTax,'unclassified')),:)=[];
        getTax(find(contains(getTax,' sp')),:)=[];
        getTax(find(contains(getTax,' bacterium')),:)=[];
        getTax(find(cellfun(@isempty,getTax)),:)=[];

        dataByTaxon(2:size(getTax)+1,1)=getTax;
        for j=1:length(getTax)
            getStrains=infoFile(find(strcmp(infoFile(:,taxCol),getTax{j})),1);
            taxonSummary=[];
            for k=1:length(getStrains)
                findModel=find(strcmp(data(:,1),getStrains{k}));
                taxonSummary(k,1:size(data,2)-1)=cell2mat(data(findModel,2:end));
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
        % then normalize the data to the highest abundance value for each
        % subsystem
        for j=2:size(dataByTaxon,2)
            maxAll=max(cell2mat(dataByTaxon(2:end,j)));
            for k=2:size(dataByTaxon,1)
                dataByTaxon{k,j}=dataByTaxon{k,j}/maxAll;
            end
        end
        for j=2:size(dataByTaxon,1)
            for k=2:size(dataByTaxon,2)
                if isnan(dataByTaxon{j,k})
                dataByTaxon{j,k}=0;
                end
            end
        end
        % limit to the 50 most abundant subsystems
        rowlabels=dataByTaxon(:,1);
        dataByTaxon(:,1)=[];
        cSums = sum(cell2mat(dataByTaxon(2:end,2:end)));
        [A,B] = sort(cSums,'descend');
        dataByTaxon=horzcat(rowlabels,dataByTaxon(:,B(1:50)));

        % plot the data
        cgo = clustergram(cell2mat(dataByTaxon(2:end,2:end))',...
            'RowLabels', dataByTaxon(1,2:end),...
            'ColumnLabels', dataByTaxon(2:end,1),...
            'ColumnLabelsRotate', 45, ...
            'Cluster', 'all', ...
            'symmetric','False', ...
            'colormap', 'jet' ...
            );
        h = plot(cgo);
        set(h,'TickLabelInterpreter','none');
        colorbar(h)
        print(['SubsystemsByTaxon' filesep paths{p,1} '_' taxLevels{i}],'-dpng','-r300')
        writetable(cell2table(dataByTaxon),['SubsystemsByTaxon' filesep paths{p,1} '_' taxLevels{i}],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

        % export taxonomic information file for R plot
        if i<6
            TaxonomyReduced=infoFile;
            taxonCol = find(strcmp(TaxonomyReduced(1, :), taxLevels{i}));
            TaxonomyReduced(:,1:taxonCol-1)=[];
            % remove duplicate entries
            [C,IA] = unique(TaxonomyReduced(:,1),'stable');
            TaxonomyReduced=TaxonomyReduced(IA,:);
            TaxonomyReduced(find(strcmp(TaxonomyReduced(:,1),'N/A')),:)=[];
            TaxonomyReduced(find(contains(TaxonomyReduced(:,1),'unclassified')),:)=[];
            TaxonomyReduced(find(contains(TaxonomyReduced(:,1),' bacterium')),:)=[];
            TaxonomyReduced(find(cellfun(@isempty,TaxonomyReduced(:,1))),:)=[];
            writetable(cell2table(TaxonomyReduced),['SubsystemsByTaxon' filesep paths{p,1} '_' taxLevels{i} '_taxonomy'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
        end
    end
end
