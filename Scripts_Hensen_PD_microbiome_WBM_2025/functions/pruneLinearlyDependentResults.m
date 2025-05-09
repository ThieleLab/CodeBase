function [prunedFluxes, linearDepMetaboliteTable, prunedFluxPath, linDepMetTablePath] = pruneLinearlyDependentResults(fluxPath, metLocationInfoFilePath, savedir)
% PRUNELINEARLYDEPENDENTRESULTS Removes linearly dependent metabolic fluxes
%   This function analyzes metabolic flux data to identify and remove linearly
%   dependent reactions, considering their presence in microbiome compartments.
%
%   NOTE! This function explicitly considers the microbiome-scaled flux
%   results with removed reactions due to low sample counts.
%
% Inputs:
%   fluxPath - Path to CSV file containing flux data
%   metLocationInfoFilePath - Path to CSV file with metabolite location information
%   savedir - Directory to save output files
%
% Outputs:
%   prunedFluxes - Table with unique (non-linearly dependent) flux results
%   linearDepMetaboliteTable - Table showing linear dependency relationships
%   prunedFluxPath - Path to saved pruned fluxes CSV file
%   linDepMetTablePath - Path to saved linear dependency table CSV file

% Input validation
validateattributes(fluxPath, {'char', 'string'}, {'nonempty'}, mfilename, 'fluxPath');
validateattributes(metLocationInfoFilePath, {'char', 'string'}, {'nonempty'}, mfilename, 'metLocationInfoFilePath');
validateattributes(savedir, {'char', 'string'}, {'nonempty'}, mfilename, 'savedir');

% Ensure save directory exists
if ~exist(savedir, 'dir')
    mkdir(savedir);
end

% Load and preprocess flux data. Table 7 contains the microbiome-scaled
% results.
fluxesOri = readtable(fluxPath, 'VariableNamingRule', 'preserve');

% Remove row column if present
if any(matches(fluxesOri.Properties.VariableNames,'sex','IgnoreCase',true))
    fluxesOri.Row = [];
end

% Remove carnitione and stearylcarnitine from the results
fluxesOri.("DM_stcrn[bc]") = [];
fluxesOri.("DM_crn[bc]") = [];

% Remove ID and sex columns
metaColumns = {'ID', 'Sex'};
fluxes = removevars(fluxesOri, metaColumns);
rxnNames = fluxes.Properties.VariableNames;

% Create flux matrix
fluxMatrix = table2array(fluxes);

% Set nan to zero to ensure that linear dependence can be calculated.
%fluxMatrix(isnan(fluxMatrix))=0;

% This for loop quantifies linear dependence by calculating the euclidean
% norm x of the difference between v1 and v2 times lambda, the ratio between 
% v2 and v1, i.e., x = ||v1 - lambda*v2||. Note that v1 and v2 are
% linearly dependent if ||v1 - lambda*v2|| = 0 

% Preallocate table for quantifying linear dependence
linDepTable = array2table(nan(size(fluxMatrix,2),size(fluxMatrix,2)),'RowNames',rxnNames','VariableNames',rxnNames);

for i=1:size(fluxMatrix,2)
    for j = 1:size(fluxMatrix,2)
    
    if i ~= j
        % Set up function input:
        vector1 = fluxMatrix(:,i);
        vector2 = fluxMatrix(:,j);

        % Find rows where either vector has a NaN value
        rowsToRemove = isnan(vector1) | isnan(vector2);
        
        % Remove rows with NaN values from both vectors
        vector1(rowsToRemove) = [];
        vector2(rowsToRemove) = [];

        % normalise the vectors to make them robust to scale differences. 
        norm_vector1 = vector1 / norm(vector1);
        norm_vector2 = vector2 / norm(vector2);
        
        % Calculate lambda by solving a linear system of equations A*lambda = B
        lambda = norm_vector2 \ norm_vector1; 
        
        % Calculate the euclean norm
        linNorm = norm(norm_vector1 - lambda * norm_vector2);
        
        % Add linear norm to table
        linDepTable{i,j} = linNorm;
    end
    
    end
end

% Apply threshold to identify linearly dependent fluxes
adjMatrix = table2array(linDepTable);

% Find identical distributions
tol = 1e-6; % Set tolerance to one in a million
adjMatrix(adjMatrix<tol) = 1;
adjMatrix(adjMatrix~=1) = 0;

% Filter on linearly dependent metabolites
colsAndRowsToRemove = all(adjMatrix==0); % Note that the columns and rows are identical

% Find linearly dependent networks
linDepTable1 = array2table(adjMatrix,'RowNames',rxnNames,'VariableNames',rxnNames);
linDepTable1(colsAndRowsToRemove',:) = [];
linDepTable1(:,colsAndRowsToRemove) = [];

% Create graph and analyze connectivity
G = graph(table2array(linDepTable1), linDepTable1.Properties.VariableNames);
nodeDegrees = degree(G);
figure;
p = plot(G);
p.Interpreter = 'none';

% Find connected components (subnetworks)
subnetworks = conncomp(G);

% Create subnetwork list with relevant information
subnetwork_list = table(subnetworks', G.Nodes.Name, nodeDegrees, ...
    'VariableNames', {'Subnetwork', 'Reaction', 'Degree'});
subnetwork_list = subnetwork_list(subnetwork_list.Degree ~= 0, :);
subnetwork_list.Degree = [];


% Process metabolite location information
metaboliteLocation = readtable(metLocationInfoFilePath,'VariableNamingRule','preserve');
metaboliteLocation.Metabolite = append('DM_', metaboliteLocation.Metabolite, '[bc]');

% Add microbiome presence information
[~, ~, ib] = intersect(subnetwork_list.Reaction, metaboliteLocation.Metabolite, 'stable');
subnetwork_list.Microbiome = metaboliteLocation.Microbiome(ib);

% If one or more metabolites are present in the microbiome compartment AND
% if one or more metabolites are NOT present in the microbiome compartment,
% REMOVE the metabolites that are NOT present in the microbiome
% compartment.

% Now for each group, identify metabolites to keep. 
[groups,groupNames] = findgroups(subnetwork_list.Subnetwork);

subnetwork_list1 = subnetwork_list;
subnetwork_list1.Keep = zeros(size(subnetwork_list1,1),1);

for i=1:length(groupNames)
    groupIndex = groups == groupNames(i);
    
    % Find metabolites in group
    groupInfo = subnetwork_list(groupIndex, :);
    
    if all(groupInfo.Microbiome,true) % all metabolites in the microbiome
        subnetwork_list1.Keep(groupIndex) = true(length(groupInfo.Microbiome),1);
    end
    
    if ~any(groupInfo.Microbiome,true) % none of the metabolites in the microbiome
        subnetwork_list1.Keep(groupIndex) = true(length(groupInfo.Microbiome),1);
    end
    
    if sum(groupInfo.Microbiome) > 0 && ... % one or more metabolites in the microbiome AND
       sum(groupInfo.Microbiome) < length(groupInfo.Microbiome) % one or more metabolites not in the microbiome 
            subnetwork_list1.Keep(groupIndex) = groupInfo.Microbiome; % Identify metabolites to keep
    end
end

% Create output table
linearDepMetaboliteTable = sortrows(subnetwork_list1,'Subnetwork');

% After removing all lindep metabolites that were not in the microbiota
% lumen, we will remove metabolites based on which metabolite in the
% subnetwork is most downstream. If that is hard to say, we will remove the
% metabolites in a subnetwork with the lowest number of nan values. 

% First, we will add the number of non-nan values for each metabolite
naVals = varfun(@(x) sum(isnan(x)), fluxes);
naVals.Properties.VariableNames = erase(naVals.Properties.VariableNames,'Fun_');
naVals = stack(naVals,1:width(naVals),'NewDataVariableName','nanVals','IndexVariableName','Reaction');
naVals.Reaction = cellstr(naVals.Reaction);
linearDepMetaboliteTable = outerjoin(linearDepMetaboliteTable,naVals,'MergeKeys',true,'Type','left');


% I can actually use the clusters to reconstruct the union of flux results
% for the entire cluster. This can be done by 1) finding the mapping
% function for each reaction pair in the cluster and normalising the fluxes
% according the to the vector with the lowest values. 2) Generate a new
% flux vector that imputes the missing values with the scaled values of
% linearly dependent reactions. 

% Find the subnetworks in the fluxes table
linDepFluxes = struct;
for i=1:max(linearDepMetaboliteTable.Subnetwork)
    rxnsToCheck = linearDepMetaboliteTable.Reaction(linearDepMetaboliteTable.Subnetwork==i);
    linDepFluxes.(strcat("subnet_",string(i))) = [fluxesOri(:,metaColumns) fluxes(:,rxnsToCheck')];
end

% Now, I want to scale all values in the linearly dependent matrices.
% First, I want to find the scaling factor for each vector and apply that.

numSubnets = max(linearDepMetaboliteTable.Subnetwork);
constructedClusterFluxes = table('Size', ...
    [height(fluxes),2+numSubnets], ...
    'VariableTypes',[{'string','string'} repmat({'double'},1,numSubnets)], ...
    'VariableNames',[{'ID'},{'Sex'} cellstr(arrayfun(@(i) append("subset_",string(i)),1:numSubnets))]);
constructedClusterFluxes.ID = fluxesOri.ID;
constructedClusterFluxes.Sex = fluxesOri.Sex;

% Add custom table property that includes the associated reactions
constructedClusterFluxes = addprop(constructedClusterFluxes,{'Reactions'},{'variable'}); 
constructedClusterFluxes.Properties.CustomProperties.Reactions = cell(width(constructedClusterFluxes),1);

for i=1:numel(fieldnames(linDepFluxes))

    % Take a subset of linearly dependent fluxes and create a matrix
    linDepTest = linDepFluxes.(strcat("subnet_",string(i)));
    linDepMatrix = table2array(linDepTest(:,3:end));
    
    % Remove rows with NaN values in the matrix
    linDepMatrixNoNaN = linDepMatrix(all(~isnan(linDepMatrix),2),:);
    
    % Find the column with the smallest flux values. Use this column as a base.
    [~, min_index] = min(median(linDepMatrixNoNaN));
    min_col = linDepMatrixNoNaN(:,min_index);
    
    % Obtain the scaling factors compared to the column with the smallest
    % values. Note that this only works because the reactions have linearly
    % dependent fluxes.
    scalingFactors = arrayfun(@(x) linDepMatrixNoNaN(:,x) \ min_col, 1:width(linDepMatrixNoNaN));
    
    % Apply the scaling factors to each column. This operation will scale the
    % values of each reaction to that of the reaction with the smallest values.
    scaled_linDepMatrix = linDepMatrix .* scalingFactors;
    
    % Now, we can perform imputation of all NaN values in the matrix column
    % with the smallest values
    
    % Create new vector for integrated reaction fluxes with imputed missing
    % values.
    clusterCol = linDepMatrix(:,min_index);
    
    % Find the rows with nan values
    naRows = find(isnan(clusterCol));
    
    % Check if any of the other column contain a non-nan values. If yes, then
    % impute the scaled values from the other columns in clusterCol.
    for j=1:length(naRows)
        nonNanCol = find(~isnan(scaled_linDepMatrix(naRows(j),:)));
        if ~isempty(nonNanCol)
            clusterCol(naRows(j)) = scaled_linDepMatrix(naRows(j),nonNanCol(1));
        end
    end
    
    % Add the imputed vector to the linDepFluxes variable
    linDepFluxes.(strcat("subnet_",string(i))).Imputed = clusterCol;

    % Also add the imputed vector to a table with the constructed cluster
    % fluxes.
    constructedClusterFluxes.(strcat("subset_",string(i))) = clusterCol;

    % Store the associated reactions in each column
    constructedClusterFluxes.Properties.CustomProperties.Reactions{i+2} = linDepTest.Properties.VariableNames(3:end);
end

% Now, we can remove all linearly dependent metabolites and replace them by
% the constructed clusters

% Create pruned fluxes table
prunedFluxes = [fluxesOri(:, metaColumns), removevars(fluxes, subnetwork_list1.Reaction)];

% Add back the constructed fluxes
prunedFluxes = outerjoin(prunedFluxes,constructedClusterFluxes,'MergeKeys',true);


% Crucially, I will also need to do that for the biomass shadow prices.
% Next, I will need to check if the shadow price values are also linearly
% dependent to each other. If so, I can perform the same imputation
% technique. If not, I can not integrate the linearly dependent fluxes. 

% Test: Load the shadow price results of cluster 4 (The largest cluster) 
% spFolder = 'C:\Users\mspg\Documents\parkinson_recreated\resultStatistics\biomass_shadow_prices';
% spRxnNames = {dir(spFolder).name}';
% spRxnNames(1:2) = [];
% spPaths = string(append(spFolder,filesep,spRxnNames));
% spRxnNames = string(erase(spRxnNames,'.csv'));
% 
% % Find all reactions in cluster 4
% 
% i=2;
% rxnsToLoad = string(linearDepMetaboliteTable.Reaction(linearDepMetaboliteTable.Subnetwork==i));
% pathsToLoad = spPaths(matches(spRxnNames,rxnsToLoad));
% rxnsForStruct = append('X',erase(rxnsToLoad,{'DM_','[bc]'}));
% spStruct = struct;
% for j = 1:length(pathsToLoad)
%     spTable = readtable(pathsToLoad(j),'VariableNamingRule','preserve');
%     spTable = renamevars(spTable,'Row','ID');
%     %spTable(:,2:end) = fillmissing(spTable(:,2:end),'constant',0);
%     spStruct.(rxnsForStruct(j)) = spTable;
% end

% If the fluxes are linearly dependent, the microbiome shadow prices are
% also linearly dependent between microbial species and samples. The nan
% values can be imputed for the shadow prices.


% Save results
prunedFluxPath = fullfile(savedir, 'uniqueFluxResults.csv');
linDepMetTablePath = fullfile(savedir, 'linearDepMetabolites.csv');

writetable(prunedFluxes, prunedFluxPath);
writetable(linearDepMetaboliteTable, linDepMetTablePath);
end
