function [corrTable, corrTableTopMicrobes, filePaths] = correlateFluxWithRa(fluxPath, fluxLimTablePath, mappedMicrobePath, saveDir)
% Correlate relative abundances with the predicted fluxes

% Load the processed fluxes and process sample IDs
fluxes = readtable(fluxPath,'VariableNamingRule','preserve');
%fluxes.ID = erase(fluxes.ID,{'muWBM_','_male','_female'});
fluxes.Sex = [];

% Load the fluxLimiter table
fluxLimTable = readtable(fluxLimTablePath,'VariableNamingRule','preserve','ReadRowNames',true);
fluxLimTable.Row = erase(fluxLimTable.Row,'pan');

% Filter the fluxes table on the reactions in fluxLimTable
rxnsOfInterest = fluxLimTable.Properties.VariableNames;
metabolites = renameVmhToMetName(rxnsOfInterest, true);

fluxes = fluxes(:,[{'ID'}, metabolites]);
fluxes.Properties.VariableNames = [{'ID'}, rxnsOfInterest];

% Load the relative abundances
microbiome = readtable(mappedMicrobePath,'VariableNamingRule','preserve','ReadRowNames',false);
microbiome = renamevars(microbiome,'Row','ID');    
microbiome = removevars(microbiome,'Sum of taxa');

% Set all nan values to zero
microbiome(:,2:end) = fillmissing(microbiome(:,2:end),'constant',0);

% Ensure that the samples are in the correct order between the fluxes and
% relative abundances
[~,ia,ib] = intersect(fluxes.ID,microbiome.ID,'stable');
fluxes = fluxes(ia,:);
microbiome = microbiome(ib,:);

% Correlate flux values with relative abundances
[RHO, ~] = corr(fluxes{:,2:end}, microbiome{:,2:end}, 'type', 'Spearman','Rows','pairwise');

% transform matrix to correlation table
corrTable = array2table(RHO',"RowNames",microbiome.Properties.VariableNames(2:end),'VariableNames',fluxes.Properties.VariableNames(2:end));

% Remove correlation coefficients of microbial species that could not
% influence its associated metabolite

% Get the associated microbial species for each reaction
fluxMicrobes = cellfun(@(x) fluxLimTable.Row(fluxLimTable.(x)==1),rxnsOfInterest,'UniformOutput',false);

% Set the correlation coefficients for each microbe to zero if not
% associated with the microbial species.
for i = 1:length(rxnsOfInterest)
    corrTable.(rxnsOfInterest{i})(~ismember(corrTable.Row,fluxMicrobes{i})) = 0;
end

% Filter on microbial species that were not analysed
microbesToKeep = any( abs(table2array(corrTable))>0 ,2);
corrTableTopMicrobes = corrTable(microbesToKeep,:);

% Set all zeros to nans
corrMatrix = table2array(corrTableTopMicrobes);
corrMatrix(corrMatrix==0) = nan;
corrTableTopMicrobes{:,:} = corrMatrix;

% Remove

% Find the top absolute correlation coefficients and remove all microbial
% species that are not in the top for any of the metabolites.

if 0
    % Find the indices of microbial species to keep
    [~,microbesToKeep] = cellfun(@(x) maxk(abs(corrTable.(x)), top),rxnsOfInterest,'UniformOutput',false);
    
    % Prune table
    corrTableTopMicrobes = corrTable(unique(vertcat(microbesToKeep{:}),'stable'),:);
end

% Sum the microbial relative abundances for each set of microbial species
summedRA = cellfun(@(x) sum(microbiome{:,x},2),fluxMicrobes,'UniformOutput',false);

% Create matrix of summed relative abundances
summedRA = [summedRA{:}];

% Correlate microbe sets with the fluxes
fluxVals = fluxes{:,2:end};
RHO = arrayfun(@(x) corr(fluxVals(:,x), summedRA(:,x), 'type', 'Spearman','Rows','pairwise'), 1:size(fluxVals,2),'UniformOutput',true);
corrAllMicrobes = array2table(RHO,"RowNames",{'Flux-associated taxa'},'VariableNames',fluxes.Properties.VariableNames(2:end));

% Put the flux-associated microbe set to the last row
corrTable = [corrTable ; corrAllMicrobes];

% Save results
filePaths.fluxCorrPath = fullfile(saveDir,'fluxRACorr.csv');
writetable(corrTable,filePaths.fluxCorrPath,'WriteRowNames',true);

% Save pruned table
corrTableTopMicrobes.Properties.VariableNames = renameVmhToMetName(corrTableTopMicrobes.Properties.VariableNames);
filePaths.prunedFluxCorrPath = fullfile(saveDir,'fluxRACorrTopMicrobes.csv');
writetable(corrTableTopMicrobes,filePaths.prunedFluxCorrPath,'WriteRowNames',true);
end