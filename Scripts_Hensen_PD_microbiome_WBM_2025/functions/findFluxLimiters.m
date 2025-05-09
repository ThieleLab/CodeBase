function [fmAssociations, fmLinksSummary, filePathIntersection, filePathStats] = findFluxLimiters(shadowPriceDir, rxnsOfInterest, saveDir, microbesOfInterest, metadataPath)
% Identify flux-limiting microbial species for metabolites of interest based on 
% shadow price data and selected metabolites. This function processes flux shadow 
% price data and calculates flux  limiting statistics for each metabolite.

if nargin<4
    microbesOfInterest = table();
end

if nargin<5
    samplesToFilter = {};
else
    samplesToFilter = readcell(metadataPath,'Range','A2:A9999'); % Only read the same IDs (row 2-end, column 1)
end

% Find paths to files with shadow price values
fileNames = {dir(shadowPriceDir).name};
fileNames(1:2) = [];

% Filter on metabolites of interest
fileNamesFiltered = fileNames(contains(fileNames,rxnsOfInterest))';

% Get reactions of interest
reactionsOfInterest = erase(fileNamesFiltered,'.csv');

% Sort files according to the rxnsOfInterest
[~,~,ib] = intersect(rxnsOfInterest,reactionsOfInterest,'stable');
fileNamesFiltered = fileNamesFiltered(ib,:);
reactionsOfInterest = reactionsOfInterest(ib);

% Find the full paths for these files
fluxLimFiles = string(append(shadowPriceDir, filesep, fileNamesFiltered));

% Process vmh ids to metabolite names
metabolites = renameVmhToMetName(reactionsOfInterest);

% Load shadow prices for reactions of interest and remove pan indications
rFiles = @(x) readtable(x,'VariableNamingRule','preserve');
rmPan = @(x) renamevars(x,x.Properties.VariableNames,erase(x.Properties.VariableNames,'pan'));
applyFun = @(x) rmPan(rFiles(x));
rxnBioShadowPrices = arrayfun(@(x) applyFun(x),fluxLimFiles,'UniformOutput',false);

% Remove the last column if present
if ismember({'Sum of taxa'},rxnBioShadowPrices{1}.Properties.VariableNames)
    rxnBioShadowPrices = cellfun(@(x) removevars(x,'Sum of taxa'), rxnBioShadowPrices,'UniformOutput',false);
end

% Remove samples not in the metadata variable
if ~isempty(samplesToFilter)
    rxnBioShadowPrices = cellfun(@(x) x(contains(x{:,1},samplesToFilter), :), rxnBioShadowPrices,'UniformOutput',false);
end

if ~isempty(microbesOfInterest)
    disp('Get summary statistics on microbial contributors for user-defined metabolite list.')
    % Get metabolites to filter
    metsToFilter = reshape(erase(microbesOfInterest.Properties.VariableNames,'pan'),[],1); % Ensure that the metabolites are a vertical array
    
    % Sort metsToFilter according to the previously found metabolites
    [~,~,ib] = intersect(metabolites,metsToFilter,'stable');
    metsToFilterSorted = metsToFilter(ib);
    
    % Get microbes to filter
    microbesToFilter = cellfun(@(x) microbesOfInterest.Row(microbesOfInterest.(x)==1), metsToFilterSorted,'UniformOutput',false);
    
    % Filter the shadow price tables on the microbes of interest
    rxnBioShadowPrices = cellfun(@(x,y) x(:,[{'Row'} y']), rxnBioShadowPrices,microbesToFilter,'UniformOutput',false);
end

% Remove ID information and obtain the value matrices
rxnBioSParray = cellfun(@(x) table2array(x(:,2:end)),rxnBioShadowPrices,'UniformOutput',false);

% Find the microbial species with non-zero shadow prices 
rxnBioSParrayLogical = cellfun(@(x) x~=0 & ~isnan(x), rxnBioSParray,'UniformOutput',false);

% Calculate the total number of microbial flux limiters and summary
% statistics for the mean and SD.
fun = @(x) [mean(sum(x,2)),std(sum(x,2)) sum(any(x))];
summaryMatrix = cell2mat(vertcat(cellfun(fun, rxnBioSParrayLogical,'UniformOutput',false)));

% Create table
fmLinksSummary = array2table(summaryMatrix,'VariableNames',{'Mean microbe count','SD microbe count','Total microbes'},'RowNames',metabolites);

% Now, lets create an intersection table for each metabolite from rxnBioSParrayLogical
if isempty(microbesOfInterest)
    microbeMetaboliteMatrix = cell2mat(cellfun(@(x) any(x),rxnBioSParrayLogical,'UniformOutput',false))';

    % Create table
    microbes = rxnBioShadowPrices{1}.Properties.VariableNames(2:end)';
    fmAssociations = array2table(microbeMetaboliteMatrix,'RowNames',microbes,'VariableNames',metabolites);

else
    fmAssociations = microbesOfInterest;
    filePathIntersection = '';
end


if isempty(microbesOfInterest)
    % Save summary stats
    filePathStats = fullfile(saveDir,'fluxLimiterStats.csv');
    writetable(fmLinksSummary, filePathStats,'WriteRowNames',true);
    
    % Save intersection table
    filePathIntersection = fullfile(saveDir,'fluxLimiterTable.csv');
    writetable(fmAssociations,filePathIntersection,'WriteRowNames',true);
else
    % Save summary stats
    filePathStats = fullfile(saveDir,'fluxLimitersFilteredStats.csv');
    writetable(fmLinksSummary, filePathStats,'WriteRowNames',true);
end

end