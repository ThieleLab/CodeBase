function [reorderedTable, meanJDist] = analyseSubsetSimilarities(taxonomyTable, plotFig,intersectionTable)
%% Calculate jaccard similarities


if nargin<3
    intersectionTable = '';
end

if ~isempty(intersectionTable)
    metabolites = intersectionTable.Properties.VariableNames;    
    binaryMatrix = table2array(intersectionTable);
else
    data = cell(height(taxonomyTable),2);
    data(:,1) = cellstr(taxonomyTable.Taxon);
    data(:,2) = cellstr(taxonomyTable.Metabolite);
    
    % Extract unique taxa and metabolites
    taxa = unique(data(:, 1));
    metabolites = unique(data(:, 2));
    
    % Create a binary matrix representing Taxon presence for each metabolite
    binaryMatrix = zeros(length(metabolites), length(taxa));
    for i = 1:length(metabolites)
        metabolitetaxa = data(strcmp(data(:, 2), metabolites{i}), 1);
        for j = 1:length(taxa)
            binaryMatrix(i, j) = any(strcmp(metabolitetaxa, taxa{j}));
        end
    end
end


% Compute pairwise Jaccard similarities between metabolites
numMetabolites = length(metabolites);
jaccardIndices = zeros(numMetabolites, numMetabolites);
for i = 1:numMetabolites
    for j = 1:numMetabolites
        intersection = sum(binaryMatrix(i, :) & binaryMatrix(j, :));
        union = sum(binaryMatrix(i, :) | binaryMatrix(j, :));
        jaccardIndices(i, j) = (intersection / union);
    end
end

% Perform hierarchical clustering based on the Jaccard similarities
linkages = linkage(squareform(1-jaccardIndices), 'average');


% Reorder the metabolites based on the clustering
dendroOrder = optimalleaforder(linkages, 1-jaccardIndices);

% Reorder the resultTable and similarities matrix
reorderedsimilarities = jaccardIndices(dendroOrder, dendroOrder);
reorderedMetabolites = metabolites(dendroOrder);

% Create a new table with the reordered metabolites
reorderedTable = array2table(reorderedsimilarities, 'VariableNames', reorderedMetabolites, 'RowNames', reorderedMetabolites);

% Find the average dissimilarity for each metabolite
j1 = reorderedsimilarities;
j1(eye(size(j1))==1) = nan;
meanJDist = table(reorderedTable.Properties.RowNames,mean(j1,'omitnan')',std(j1,'omitnan')');
meanJDist.Properties.VariableNames = {'Metabolite','Mean Jaccard similarity','SD jaccard similarity'};

if plotFig
    figure;
    heatmap(reorderedMetabolites, reorderedMetabolites, reorderedsimilarities, 'Colormap', parula, 'ColorLimits', [0 1]);
    jaccDist = reshape(triu(jaccardIndices),[],1);
    jaccDist(jaccDist==0)=[];
    sub = append("Mean = ",string(round(mean(jaccDist,'omitnan'),3)),", SD = ", string(round(std(jaccDist,'omitnan'),3)));
    title({'Heatmap of Clustered Jaccard similarities for Metabolites', char(sub)});
    xlabel('');
    ylabel('');
end