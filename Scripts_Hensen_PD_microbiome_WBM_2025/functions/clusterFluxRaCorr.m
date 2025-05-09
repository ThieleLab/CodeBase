function [bestCorrPerSizeList, minbestCorrPerSizeList, rxnsOfInterest] = clusterFluxRaCorr(paths, top, onlyCorrTopMicrobes)
%function [bestCorrPerSizeList, minbestCorrPerSizeList, preparedScatterPlotData, clusterCorr, clusterAbundancesCell, rxnsOfInterest] = clusterFluxRaCorr(paths, top, cutoff, onlyCorrTopMicrobes)

% INPUTS:
fluxPath = paths.fluxPath;
fluxMicrobeCorrPath = paths.fluxRaCorrPath;
mappedMicrobePath = paths.mappedMicrobePath;


% Read correlation data
corrTableOri = readtable(fluxMicrobeCorrPath, 'PreserveVariableNames', true);
corrTableOri.Row = strrep(corrTableOri.Row, ' ', '_');
corrTable = corrTableOri;
corrTable(matches(corrTable.Row, 'Flux-associated_taxa'),:)= [];

% Get reactions of interest
rxnsOfInterest = corrTable.Properties.VariableNames(2:end);

if onlyCorrTopMicrobes
    % Identify the top X microbial species for each metabolite
    [~,microbesToKeep] = cellfun(@(x) maxk(abs(corrTable.(x)), top),rxnsOfInterest,'UniformOutput',false); % Get row indices
    %[~,microbesToKeep] = cellfun(@(x) mink(corrTable.(x), top),rxnsOfInterest,'UniformOutput',false); % Get row indices
    topMicrobeLists = cellfun(@(x) corrTable.Row(x), microbesToKeep,'UniformOutput',false); % Get taxon names

else 
    % Test all flux-associated microbes
    microbesToKeep = cellfun(@(x) find(corrTable.(x)~=0), rxnsOfInterest,'UniformOutput',false); % Get row indices
    topMicrobeLists = cellfun(@(x) corrTable.Row(x), microbesToKeep,'UniformOutput',false); % Get taxon names
end

% Load relative abundance information
relAbun = readtable(mappedMicrobePath, 'PreserveVariableNames', true);
relAbun.Properties.VariableNames = strrep(relAbun.Properties.VariableNames, 'pan', '');

% Replace nan values with zeros
relAbun{:,2:end} = fillmissing(table2array(relAbun(:,2:end)), 'constant', 0);

% Finally, load flux data into memory
fluxValues = readtable(fluxPath,'VariableNamingRule','preserve');

% Rename VMH IDs of metabolites of interest to metabolite names
fluxValues.Properties.VariableNames = renameVmhToMetName(fluxValues.Properties.VariableNames);
%%
bestCorrPerSizeList = cell(1,length(rxnsOfInterest));
minbestCorrPerSizeList = cell(1,length(rxnsOfInterest));

for i=1:length(rxnsOfInterest)
    topMicrobes = topMicrobeLists{i};
    rxnName = rxnsOfInterest{i};

    [clusterNames, clusterSize, clusterAbundances] = createClusters(top, topMicrobes, relAbun);

    rxnClusterCorr = correlateClusters(fluxValues, clusterAbundances, rxnName, clusterNames, clusterSize);
    clearvars clusterNames clusterSize clusterAbundances
    
    
    % Get the top X clusters
    [bestCorrPerSizeList{i},minbestCorrPerSizeList{i}] = extractTopClustersPerSize(rxnClusterCorr,rxnName);
    clearvars clusterCorr
end
% %%
% % Create clusters of summed relative abundances
% [clusterNamesCell, clusterSizeCell, clusterAbundancesCell] = cellfun(@(topMicrobes) createClusters(top, topMicrobes, relAbun), topMicrobeLists,'UniformOutput',false);
% %%
% 
% % Correlate fluxes with microbe clusters
% clusterCorr = cellfun(@(clusterAbundances, rxnName, clusterNames, clusterSize) correlateClusters(fluxValues, clusterAbundances, rxnName, clusterNames, clusterSize), ...
%     clusterAbundancesCell, rxnsOfInterest, clusterNamesCell, clusterSizeCell ,'UniformOutput',false);
% 
% %%
% % Get the top X clusters
% [bestCorrPerSizeList,minbestCorrPerSizeList] = cellfun(@(rxnClusterCorr,rxnName) extractTopClustersPerSize(rxnClusterCorr,rxnName),... 
% clusterCorr, rxnsOfInterest,'UniformOutput',false);
% %%
% % Extract the flux values and corresponding top combined microbes
% preparedScatterPlotData = cellfun(@(bestCorrResults,clusterAbundances,rxnName) prepareTopClustForScatterPlot(fluxValues, bestCorrResults, clusterAbundances, rxnName, cutoff),... 
% bestCorrPerSizeList, clusterAbundancesCell, rxnsOfInterest,'UniformOutput',false);

end

function [clusterNames, clusterSize, clusterAbundances] = createClusters(top, topMicrobes, relAbun)
% Pre-calculate number of clusters for each size
clusterCounts = arrayfun(@(x) nchoosek(length(topMicrobes),x), 1:top);
numClusters = sum(clusterCounts);

% Pre-allocate arrays
clusterNames = strings(numClusters, top);
clusterAbundances = zeros(height(relAbun), numClusters, 'single');
clusterSize = zeros(numClusters, 1);

% Convert relAbun to numeric array once (outside loops)
numericAbun = single(table2array(relAbun(:, topMicrobes)));

% Initialize tracking variable
currentRow = 1;

% Generate combinations using a more efficient approach
for k = 1:top
    % Pre-allocate temporary arrays for batch processing
    combIdx = nchoosek(1:length(topMicrobes), k);
    batchSize = size(combIdx, 1);
    rowIndices = currentRow:(currentRow + batchSize - 1);
    
    % Set cluster sizes in batch
    clusterSize(rowIndices) = k;
    
    % Process names in batch of size k
    for j = 1:k
        clusterNames(rowIndices, j) = topMicrobes(combIdx(:, j));
    end

    % Create a sparse logical matrix for this batch
    batchIndicator = false(length(topMicrobes), batchSize);
    for i = 1:batchSize
        batchIndicator(combIdx(i,:), i) = true;
    end
    
    % Calculate abundances using matrix multiplication.
    % Note that for every matrix multiplication c(i,j) = a(i,k)*b(k,j) is defined as 
    % for 1:k -> c(i,j) = sum( a(i,k) * b(k,j) ). This multiplication
    % speeds up the summation of the microbial abundances.
    clusterAbundances(:, rowIndices) = numericAbun * batchIndicator;
    
    % Update current row position
    currentRow = currentRow + batchSize;
end
end

function clusterCorr = correlateClusters(fluxValues, clusterAbundances, rxnName, clusterNames, clusterSize)
% Get reaction fluxes for metabolite
rxnFlux = fluxValues.(rxnName);
disp(rxnName)
% Remove nan values in the fluxes and abundances for pairwise correlations
rowsToRemove = ~isnan(rxnFlux);
rxnFlux = rxnFlux(rowsToRemove);
clusterAbundances = clusterAbundances(rowsToRemove,:);

% Perform pearson correlation. Note that a spearman correlation would have
% been a bit more appropriate, however the spearman algorithm is much
% (much) slower and a linear regression of spearman and pearson results found that
% between 75 and 90% of the variance in the spearman results was captured
% in the pearson results. That is good enough for me. 
RHO = corr(rxnFlux, clusterAbundances,'type','Spearman');

% Create table that combines the cluster names with the cluster sizes and
% the calculated correlations. 
clusterCorr = array2table(clusterNames);
clusterCorr.Size = clusterSize;
clusterCorr.(rxnName) = RHO';

end

function [bestCorrPerSize,minbestCorrPerSize] = extractTopClustersPerSize(rxnClusterCorr,rxnName)
% Extract the top clusters for each cluster size

% Get the absolute correlation coefficients
% absCorrCoef = abs(rxnClusterCorr.(rxnName));
absCorrCoef = double(rxnClusterCorr.(rxnName));
% Find the group sizes
groupIds = rxnClusterCorr.Size;
% Get number of group sizes
numGroups = max(groupIds);

% Create sparse matrix where each column represents a group
% and rows represent the original indices
S = sparse(1:numel(absCorrCoef), rxnClusterCorr.Size, absCorrCoef, numel(absCorrCoef), numGroups);

% Find maximum value for each group
[~, maxIndices] = max(S, [], 1);
[~, minIndices] = min(S, [], 1);

% Extract the clusters with the highest correlations for each cluster size
bestCorrPerSize = rxnClusterCorr(maxIndices',:);
minbestCorrPerSize = rxnClusterCorr(minIndices',:);

 
% Perform the following steps to prepare the correlation tables
reVarsMicrobes = @(x) renamevars(x,x.Properties.VariableNames,replace(x.Properties.VariableNames,'clusterNames','Microbe_')); % Rename "cluster" to "Microbe"
reVar = @(x) renamevars(x,{'Size',rxnName},{'Cluster_size','Spearman_RHO'}); % Rename correlation column
addMet = @(x) addvars(x,repmat(cellstr(rxnName),height(x),1),'Before',1,'NewVariableNames','Metabolite'); % Specify associate metabolite as a separate column
addImpr = @(x) addvars(x,[nan; diff(abs(x.Spearman_RHO))],'After','Spearman_RHO','NewVariableNames','Absolute_RHO_improvement'); % Add the improved rho values for each increase in size
applyF = @(x) addImpr(addMet(reVar(reVarsMicrobes(x)))); % Combine defined table transformations

% Process top correlation tables 
bestCorrPerSize = applyF(bestCorrPerSize);
minbestCorrPerSize = applyF(minbestCorrPerSize);

% Add correlation direction information
bestCorrPerSize = addvars(bestCorrPerSize,repmat('Max',height(bestCorrPerSize),1),'After','Metabolite','NewVariableNames','Direction');
minbestCorrPerSize = addvars(minbestCorrPerSize,repmat('Min',height(minbestCorrPerSize),1),'After','Metabolite','NewVariableNames','Direction');

% Add original cluster row index
bestCorrPerSize.clusterIndex = maxIndices';
minbestCorrPerSize.clusterIndex = minIndices';
end

function preparedValsTable = prepareTopClustForScatterPlot(fluxValues, bestCorrResults, clusterAbundances, rxnName, cutoff)
% Create cell array with the flux values and the selected clusters for
% visualisation

% Find the names of the last clusters that improve on the correlations by
% 5% or more. In  other words, find the last cluster size that improves 5%
% or more.
if max(bestCorrResults.Absolute_RHO_improvement)<0.05
    cutoff = max(bestCorrResults.Absolute_RHO_improvement);
end

topClustIndex = find( bestCorrResults.Absolute_RHO_improvement>cutoff ,1, "last");
clustIndex = bestCorrResults.clusterIndex(topClustIndex);
%[bestCorrResults(:,[1 end])]
% top clusters 
topClustName = append("cluster_",string(clustIndex));

% Get the values for the associated cluster
topClustVals = clusterAbundances(:,clustIndex);

% Create table with associated values
preparedValsTable = fluxValues(:,string(rxnName));
preparedValsTable.(string(topClustName)) = topClustVals;


% Log transform the fluxes and relative abundances for visualisation
% purposes
preparedValsTable.(rxnName) = log2(preparedValsTable.(rxnName));
preparedValsTable.(string(topClustName)) = log2(preparedValsTable.(string(topClustName)));

end