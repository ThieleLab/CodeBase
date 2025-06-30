function [adjMatrix,filePath] = extractTopFluxMicrobeLinks(paths)
% matrix for microbes
% INPUT:
% paths structure
% 
% output:
% adjMatrix
% filePath

% Load table
topClusterCorr = readtable(paths.topClusterCorr);

% algorithm:

% 1) Remove rows where cluster size is smaller than two or larger than 4
% 2) Remove rows where Absolute_RHO_improvement < 0.05
% 3) Find the rows for the max absolute cluster size
% 4) Generate matrix for metabolite - microbe adjacency. 


% 1) Remove rows where cluster size is smaller than two or larger than 4
rowsToKeepCSize = topClusterCorr.Cluster_size>1 & topClusterCorr.Cluster_size<5;
prunedClustSize = topClusterCorr(rowsToKeepCSize,:);

% 2) Remove rows where Absolute_RHO_improvement < 0.05
rowsToKeepRhoImpr = prunedClustSize.Absolute_RHO_improvement>0.05;
prunedClustSize = prunedClustSize(rowsToKeepRhoImpr,:);

% 3) Find the rows for the max cluster size
prunedClustSize = groupfilter(prunedClustSize, "Metabolite",@(x) abs(x) == max(abs(x)),"Spearman_RHO");

% 4) Generate matrix for metabolite - microbe adjacency. 

% 4.1) Select only the metabolite names and associated microbes
metMicrobeLinks = removevars(prunedClustSize,{'Absolute_RHO_improvement','Cluster_size','Direction','Spearman_RHO'});

% 4.2) Create adjacency matrix
metMicrobeLinks = stack(metMicrobeLinks,2:6,'NewDataVariableName','Species','IndexVariableName','varName'); % pivot to tall table
metMicrobeLinks = removevars(metMicrobeLinks,'varName'); % Remove unnessecary column
metMicrobeLinks(matches(metMicrobeLinks.Species,""),:)=[]; % Remove empty rows
metMicrobeLinks = addvars(metMicrobeLinks,ones(size(metMicrobeLinks,1),1) , 'NewVariableNames','value');
adjMatrix = unstack(metMicrobeLinks,"value",'Metabolite','VariableNamingRule','preserve'); % Get adjacency matrix
adjMatrix(:,2:end) = fillmissing(adjMatrix(:,2:end),'constant',0); % Remove nan

% Save table
filePath = fullfile(paths.microbeToFlux,'topClusterMicrobePresenceTable.csv');
writetable(adjMatrix,filePath)
end

