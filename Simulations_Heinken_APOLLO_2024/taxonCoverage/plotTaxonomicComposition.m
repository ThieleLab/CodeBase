
% create pie charts/dnot charts of taxonomic composition of
% GlobalBiome/AGORA2

infoFileGB = readInputTableForPipeline('Combined_taxonomy_info.xlsx');

infoFileAGORA2 = readInputTableForPipeline('AGORA2_infoFile.xlsx');

taxLevels = {'Species','Genus','Family','Order','Class','Phylum'};

for i=1:length(taxLevels)
    taxCol=find(strcmp(infoFileGB(1,:),taxLevels{i}));
    [taxa, ~, J] = unique(infoFileGB(2:end,taxCol));
    cnt = histc(J, 1:numel(taxa));
    cnt(find(strcmp(taxa,' ')),:)=[];
    taxa(find(strcmp(taxa,' ')),:)=[];
    cnt(find(strcmp(taxa,'')),:)=[];
    taxa(find(strcmp(taxa,'')),:)=[];
    cnt(find(contains(taxa,'unclassified')),:)=[];
    taxa(find(contains(taxa,'unclassified')),:)=[];
    cnt(find(contains(taxa,' bacterium')),:)=[];
    taxa(find(contains(taxa,' bacterium')),:)=[];
    cnt(find(cellfun(@isempty,taxa)),:)=[];
    taxa(find(cellfun(@isempty,taxa)),:)=[];

    % get the 25 most abundant taxa if there are more than 20
    if length(taxa)>20
        [A,B] = sort(cnt,'descend');
        taxa=taxa(B);
        cnt=cnt(B);
        data=zeros(20,2);
        data(1:20,2) = cnt(1:20);
        taxaNew=taxa(1:20);
        [A,B] = sort(taxaNew);
        taxaNew=taxaNew(B);
        data(:,2)=data(B,2);
        taxaNew = vertcat(taxaNew,{'Other'},{'Unclassified'});
        data(21,2) = sum(cnt(21:end));
        data(22,2) = (size(infoFileGB,1)-1)-sum(cnt);

    else
        taxaNew = vertcat(taxa,{'Unclassified'});
        data=zeros(length(taxa)+1,2);
        data(1:length(taxa),2) = cnt;
        data(length(taxa)+1,2) = (size(infoFileGB,1)-1)-sum(cnt);
    end

    taxCol=find(strcmp(infoFileAGORA2(1,:),taxLevels{i}));
    [taxa, ~, J] = unique(infoFileAGORA2(2:end,taxCol));
    cnt = histc(J, 1:numel(taxa));
    cnt(find(strcmp(taxa,' ')),:)=[];
    taxa(find(strcmp(taxa,' ')),:)=[];
    cnt(find(strcmp(taxa,'')),:)=[];
    taxa(find(strcmp(taxa,'')),:)=[];
    cnt(find(contains(taxa,'unclassified')),:)=[];
    taxa(find(contains(taxa,'unclassified')),:)=[];

    for j=1:size(data,1)-2
        findTax=find(strcmp(taxa,taxaNew{j,1}));
        if ~isempty(findTax)
            data(j,1)=cnt(findTax);
        end
    end
    if length(taxa)>20
        data(size(data,1)-1,1) = sum(cnt(21:end));
    end
    data(size(data,1),1) = (size(infoFileAGORA2,1)-1)-sum(cnt);

    %     % define random colors
    %     cols=[];
    %     for j=1:length(taxaNew)
    %      cols(j,:)=[rand rand rand];
    %     end
    % use existing colormap
    cols = hsv(length(taxaNew));

    figure
    donut(data',{},cols)
    hold on
    legend(taxaNew,'Location','EastOutside')
    title(taxLevels{i})
    set(gca, 'FontSize', 12)
    axis off

    print(['Computation_Figure_2' filesep 'TaxonComposition_' taxLevels{i}],'-dpng','-r300')
end
