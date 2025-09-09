function updated = fixMassiliensisFiles(mappedMicrobePath, fluxLimDir)

% Get paths of data to fix
% Combine shadow prices for Enorma_massiliensis and Collinsella_massiliensis

% INPUT
% mappedMicrobePath = paths.mappedMicrobePath;
% fluxLimDir = paths.fluxLimDir;

% Find paths to files with shadow price values
fileNames = {dir(fluxLimDir).name};
fileNames(1:2) = [];
pathsToMicrobes = fullfile(fluxLimDir,fileNames)';
pathsToMicrobes = [cellstr(mappedMicrobePath); pathsToMicrobes];

% Run function for all paths
updated = cellfun(@(x) fixMasiliensis(x), pathsToMicrobes,'UniformOutput',true);
end

function updated = fixMasiliensis(pathToMicrobes)
% Load relative abundance data
microbiome = readtable(pathToMicrobes,'VariableNamingRule','preserve');

% Remove pan from microbe names
microbiome.Properties.VariableNames = erase(microbiome.Properties.VariableNames,'pan');

updated = false;
% Check if Collinsella_massiliensis is still present in the RA data
if any(matches(microbiome.Properties.VariableNames,'Collinsella_massiliensis'))
    % Combine the relative abundances for Enorma and Collinsella masiliensis
    microbiome.Enorma_massiliensis = sum(microbiome{:,{'Enorma_massiliensis','Collinsella_massiliensis'}},2,'omitnan');
    % Set all zeros to nan
    microbiome.Enorma_massiliensis(microbiome.Enorma_massiliensis==0) = nan;
    % Remove the collinsella results
    microbiome.Collinsella_massiliensis = [];
    % Save updated data
    writetable(microbiome,pathToMicrobes)

    updated = true;
end

end

