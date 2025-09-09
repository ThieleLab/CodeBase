function [bootMeanTable, fileName] = fluxMicrobeSensitivityAnalysis(fluxLimTablePath, mContributionDir, saveDir, bootSamp)
% Goal: Perform a sensitivity analysis to find the most 
% important microbial contributors to the fluxes and the most important
% flux consumers. 

% Load flux-associated microbes
fluxLimTable = readtable(fluxLimTablePath, 'PreserveVariableNames', true,'ReadRowNames',true);
% fluxLimTable.Properties.RowNames = erase(fluxLimTable.Properties.RowNames,'pan');
% fluxLimTable("Sum of taxa",:) = [];

% Convert metabolite names to rxn IDs
fluxLimTable.Properties.VariableNames = renameVmhToMetName(fluxLimTable.Properties.VariableNames,true);

% Get reactions to check
rxnsOfInterest = fluxLimTable.Properties.VariableNames;

%%

% Find the bootstrapped mean averages of potential microbial contributions
% for all metabolites 

% initialise parallel pool
poolobj = gcp('nocreate');
numRxns = numel(rxnsOfInterest);
if isempty(poolobj); parpool(numRxns); end


% Perform bootstrap analysis on all rxnsOfInterest
bootMeanTables = cell(1,numRxns);
parfor i = 1:numRxns
    bootMeanTables{i} = getMicrobeSensitivity(mContributionDir, rxnsOfInterest{i}, fluxLimTable, bootSamp);
end
% tic
% % Perform function on all metabolites
% bootMeanTables1 = cellfun(@(metabolite) getMicrobeSensitivity(mContributionDir, metabolite, fluxLimTable, bootSamp),rxnsOfInterest,'UniformOutput',false);
% toc
%%
% Create tall table with 
bootMeanTable = vertcat(bootMeanTables{:});

% Rename vmh ids to metabolite names
bootMeanTable.Metabolite = renameVmhToMetName(bootMeanTable.Metabolite);

% Write table to file
fileName = fullfile(saveDir,'microbeSensitivityTable.csv');
writetable(bootMeanTable,fileName)
end
% Add the average flux contributions to fluxLimTable
% fluxLimTable1 = fluxLimTable;
% for i=1:length(rxnsOfInterest)
%     [~,ia,ib] = intersect(fluxLimTable.Row,bootMeanTables{i}.Row,'stable');
%     fluxLimTable1.(rxnsOfInterest{i})(ia) = bootMeanTables{i}.Mean(ib);
% end
% 
% % Get the top 5 microbial species
% fluxLimTable1.Mean = mean(fluxLimTable1{:,:},2);
% fluxLimTable1 = sortrows(fluxLimTable1,'Mean','descend');
% 
% % Filter on the top 5 microbial species
% n = 10;
% fluxLimTable2 = fluxLimTable1(1:n,:);
% fluxLimTable2= sortrows(fluxLimTable2,'Mean','ascend');
% fluxLimTable2.Mean = [];
% 
% % Create barplot
% 
% close all
% figure;
% barh(fluxLimTable2{:,:})
% yticklabels(fluxLimTable2.Row);
% xlabel({'Predicted mean reduction in flux upon removal of microbial species'})
% set(gca,'TickLabelInterpreter','none')
% metabolites = renameVmhToMetName(fluxLimTable2.Properties.VariableNames);
% legend(metabolites,'Location','southwest');

%ylim([0 1])

function bootMeanTable = getMicrobeSensitivity(mContributionDir, metabolite, fluxLimTable, bootSamp)
% Load flux contributions
contributions = readtable(fullfile(mContributionDir,filesep, metabolite), 'PreserveVariableNames', true,'ReadRowNames',true);

rmPan = @(x) renamevars(x,x.Properties.VariableNames,erase(x.Properties.VariableNames,'pan'));
contributions = rmPan(contributions);

% Remove pan indications
%contributions.Properties.VariableNames = erase(contributions.Properties.VariableNames,'pan');

% Get reduction in flux when removing microbe
contributions{:,:} = contributions{:,:};
%contributions = varfun(@(x) x./fluxValues.(metabolite), contributions);
%contributions.Properties.VariableNames = erase(contributions.Properties.VariableNames,'Fun_');

% Set all nan values to zeros
contributions = fillmissing(contributions, 'constant', 0);

% Remove the sum of taxa column
%contributions(:,"Sum of taxa") = [];

% Filter on microbial species that are associated with the fluxes
fluxAssocMicrobes = fluxLimTable.Row(fluxLimTable.(metabolite)==1)';
contributions = contributions(:,fluxAssocMicrobes);

% Find the average flux contributions 
contributionMatrix = table2array(contributions);

% Calculate the mean flux contributions and the associated 95% confidence
% intervals
[ci,bootstat] = arrayfun(@(x) bootci(bootSamp,@mean,contributionMatrix(:,x)), 1:size(contributions,2),'UniformOutput',false);

% Get means and CI
bootMeans = [cellfun(@mean,bootstat)',[ci{:}]'];

% Generate table
bootMeanTable = array2table(bootMeans,'VariableNames',{'Mean','2.5CI','97.5CI'});

% Add taxa
bootMeanTable.Taxa = contributions.Properties.VariableNames';

% Add associated metabolite
bootMeanTable.Metabolite = cellstr(repmat(metabolite,height(bootMeanTable),1));

% Process table
bootMeanTable = movevars(bootMeanTable,'Metabolite','before',1);
bootMeanTable = movevars(bootMeanTable,'Taxa','after','Metabolite');

end


