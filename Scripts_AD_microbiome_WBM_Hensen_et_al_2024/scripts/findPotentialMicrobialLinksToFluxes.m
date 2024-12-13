function findPotentialMicrobialLinksToFluxes(mappedSpeciesPath, wbmPath, fluxPath, savePath)
% This function identifies microbiome species associated with specific samples,
% filters models based on sample IDs, and maps metabolites to microbes 
% based on their extracellular metabolite production capabilities. The
% function produces an excel file, metaboliteMicrobeLinks.xlsx, which
% contains information on which microbes could produce which metabolites
% based on their metabolic content. The function has the following
% sections:
% 1. Load and preprocess microbiome data to identify relevant sample IDs.
% 2. Filter metabolic models based on identified sample IDs.
% 3. Extract metabolite data and establish metabolite-microbe links.
% 4. Save results in a structured output file.
%
% USAGE:
%   findPotentialMicrobialLinksToFluxes(mappedSpeciesPath, wbmPath, fluxPath)
%
% INPUTS:
% - mappedSpeciesPath: Contains microbiome data with species and sample IDs.
% - fluxPath: Contains metabolite flux data for analysis.
% - wbmPath: WBM models directory. Path containing metabolic model files.
% - savePath: Path to directory where metaboliteMicrobeLinks.xlsx will be saved 
%
% OUTPUTS:
%   None    
%
% .. Author:
%       - Tim Hensen, November 2024

% Load microbiome species data from CSV file
microbiome = readtable(mappedSpeciesPath);

% Extract sample IDs (column headers from the second column onward)
IDs = microbiome.Properties.VariableNames(2:end);

% Initialize empty arrays for storing included sample IDs and species
IDsToInclude = string();
samples = string();
includedSpecies = string();

% In the next step, the minimum number of samples is defined that still
% contain all species within the sample dataset. 

% Loop through each column (sample) to identify samples containing non-zero values
for i = 2:width(microbiome)
    % Identify species associated with the sample (non-zero values)
    new = string(microbiome{microbiome{:, i} ~= 0, 1});
    
    % Update the list of species and associated sample IDs
    samplesNew = union(samples, new);
    if length(samplesNew) > length(samples)
        samples = samplesNew;
        IDsToInclude = union(IDsToInclude, IDs(i));
    end
end

% Remove the first empty entry and clean up sample ID names
IDsToInclude = IDsToInclude(2:end);
IDsToInclude = erase(IDsToInclude, 'Sample');

% Get the list of model names and IDs from the specified directory
[modelNames, modelIDs] = findModelPaths(wbmPath, '.mat');

% Match models with the identified sample IDs
[~, idx] = intersect(modelIDs, IDsToInclude, 'stable');
modelNames = modelNames(idx);

% Get metabolites to test
% Load metabolite flux results from the CSV file
fluxes = readtable(fluxPath, 'VariableNamingRule', 'preserve');

% Extract metabolite IDs (column names starting from the third column)
metaboliteList = fluxes.Properties.VariableNames(3:end)';
metabolites.VMHID = metaboliteList;

% Create an empty table to store metabolite-microbe links
metMicrobeLinks = array2table(...
    zeros(height(samples) - 1, length(metabolites.VMHID)), ...
    'RowNames', microbiome.Species', 'VariableNames', metabolites.VMHID);

% Link microbes with metabolites
% Loop through each metabolic model
for i = 1:length(modelNames)
    % Load the metabolic model
    disp(strcat("Load model:", string(i), " out of:", string(length(modelNames))))
    model = load(modelNames(i)); % Load the model data
    model = model.(string(fieldnames(model))); % Extract the model structure

    % Loop through each metabolite
    disp('Investigate metabolites')
    for j = 1:length(metabolites.VMHID)
        % Get the current metabolite ID
        met = string(metabolites.VMHID(j));
        
        % Check if the metabolite is present in the extracellular microbiome lumen
        mets = model.mets(contains(model.mets, append('_', met, '[luM]')));
        
        % Get the species names to check
        microbesToCheck = metMicrobeLinks.Properties.RowNames;
        
        % Identify species capable of excreting the metabolite
        func = @(x) any(contains(mets, x));
        id = cellfun(func, microbesToCheck);
        metMicrobeLinks{id, j} = 1; % Mark the link in the table
    end
    disp('Added metabolites')
end

% Save the results to an Excel file
filename = strcat(savePath, filesep, 'metaboliteMicrobeLinks.xlsx');
writetable(metMicrobeLinks, filename, 'WriteRowNames', true)
end