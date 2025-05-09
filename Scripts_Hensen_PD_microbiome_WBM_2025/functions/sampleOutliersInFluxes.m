function [sampleImportanceTable, explVar] = sampleOutliersInFluxes(fluxPath, metadataPath, saveDir)

% Define variable of interest for visualisation purposes
varOfInterest = 'Case_status';

% Load metadata and obtain sample case status
metadata = readtable(metadataPath,'VariableNamingRule','preserve');

% Load flux results and metadata
fluxTable = readtable(fluxPath,'VariableNamingRule','preserve');

% Remove sex information
fluxTable.Sex = [];

% Process ID information
fluxTable.ID = erase(fluxTable.ID,{'muWBM_','_male','_female'});

% Remove samples not in metadata
[~,ia,ib] = intersect(metadata.ID,fluxTable.ID,'stable');
metadata = metadata(ia,:);
fluxTable = fluxTable(ib,:);

% Create matrix
fluxNorm = normalize(table2array(fluxTable(:,2:end)));

% Perform robust PCA
warning('off')
[~, PCscores, ~, ~, explVar] = pca(fluxNorm, 'Centered', false, 'Algorithm', 'als');
warning('on')

% Calculate score importances
scoreImportance = PCscores .* explVar';

% Calculate euclidean norm
fun = @(x) norm(scoreImportance(x,:));
scoreNorms = arrayfun(fun,1:height(scoreImportance))';

sampleImportanceTable = table(fluxTable.ID,scoreNorms,metadata.(varOfInterest),'VariableNames',{'ID','norm','Disease status'});
sampleImportanceTable = sortrows(sampleImportanceTable,'norm','descend');


%%% Visualise PCA


uniqueCategories = unique(metadata.(varOfInterest));
colorMap = lines(length(uniqueCategories));


% Create plot on the explained variances
cumExplVar = cumsum(explVar);
cumExplVar = cumExplVar(cumExplVar<90);

figure;
plot(cumExplVar);
xlabel('Principal component')
ylabel('Explained variance in %')
title('PC explained variances')

% Create 3D PCA score plot

f2 = figure;
% PD samples
idx = matches(metadata.(varOfInterest),'PD');
scatter(PCscores(idx,1), PCscores(idx,2), 'MarkerEdgeColor',colorMap(1,:));
hold on;
% Controls
idx = matches(metadata.(varOfInterest),'Control');
scatter(PCscores(idx,1), PCscores(idx,2), 'MarkerEdgeColor',colorMap(2,:));
% Sample ID labels
text(PCscores(:,1), PCscores(:,2), fluxTable.ID,'FontSize',7)

% Axis titles
xlabel(strcat("PC1 explained variance: ",string(explVar(1)),"% "))
ylabel(strcat("PC2 explained variance: ",string(explVar(2)),"% "))
title('PCA of sample flux predictions')

%plot(0,0,'o','MarkerSize',200);
%plot(0,0,'o','MarkerSize',300);

% Legend
legend({'PD','Control'})

xlim([-27.5 27.5])
ylim([-27.5 27.5])

% Save figure
filePath = fullfile(saveDir,'PCA_fluxes.png');
exportgraphics(f2, filePath, "Resolution", 300)
end