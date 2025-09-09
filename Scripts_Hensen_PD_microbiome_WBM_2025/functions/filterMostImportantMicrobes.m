function [intersectionTable, bootMeanCell, filePath,fnameFullFluxContrTable] = filterMostImportantMicrobes(bootMeanTable, cumulativeFraction, saveDir)
% Find the smallest set of microbial species that would reduce the
% microbiome contribution to the fluxes by 95%

% Create cell array with means for each metabolite. 
[groups,Metabolites] = findgroups(bootMeanTable.Metabolite);
bootMeanCell = arrayfun(@(x) bootMeanTable(groups==x,:), 1:max(groups),'UniformOutput',false);

% Then, sort microbes on their sensitivity and calculate the cumulative sum
% of microbial contribution. Then, calculate the fraction of the total
% cumulative sum for each entry. 
sortM = @(x) sortrows(x,'Mean','descend');
calcCS = @(x) addvars(x, cumsum(x.Mean), cumsum(x.Mean) / sum(x.Mean),'NewVariableNames',{'cumulative_sum','fractional_cumulative_sum'});
addThreshold = @(x) addvars(x, x.fractional_cumulative_sum<cumulativeFraction, 'NewVariableNames',{'Analysed'}); 
bootMeanCell = cellfun(@(x) addThreshold(calcCS(sortM(x))), bootMeanCell,'UniformOutput',false); % Apply nested functions

% Finally, we remove all microbial species after a cutoff of 95%. 
bootMeanCellFiltered = cellfun(@(x) x(x.Analysed, :),bootMeanCell, 'UniformOutput', false);

% After removing the microbial species, we reconstruct a new intersection
% table for the flux-associated microbial species

% Translate cell array to presence table
% Get all unique microbes
microbesCell = cellfun(@(x) x.Taxa, bootMeanCellFiltered,'UniformOutput',false);
microbes = unique(vertcat(microbesCell{:}));

% Create matrix that indicates which microbes are present in each set
matrix = zeros(numel(microbes),numel(Metabolites));
for i=1:numel(Metabolites)
    matrix(matches(microbes,microbesCell{i}),i) = 1;
end

% Convert matrix to table
intersectionTable = array2table(matrix,'RowNames',microbes,'VariableNames',Metabolites);

% Save table to file
filePath = fullfile(saveDir,'mostImportantMicrobialLimiters.csv');
writetable(intersectionTable,filePath,'WriteRowNames',true);

% Write full bootMean table to file
bootMeanTable = vertcat(bootMeanCell{:});

fnameFullFluxContrTable = fullfile(saveDir,'microbeSensitivityTable.csv');
writetable(bootMeanTable,fnameFullFluxContrTable,'WriteMode','overwrite')
end



