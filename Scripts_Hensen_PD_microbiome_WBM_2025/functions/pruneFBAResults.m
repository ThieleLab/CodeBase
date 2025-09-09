function newFolder = pruneFBAResults(fbaFolder, outputFolder, metadataPath)
% This function finds the FBA solution results present in the metadata and
% copies them to a new folder. Additionally, the gf results are also
% copied.
%
% INPUT
% fbaFolder            Character array with path to FBA solutions.
%
% OUTPUT
% prunedFBApaths     Path to slimmed down FBA results
%
% AUTHOR: Tim Hensen, October 2024

% Load metadata
metadata = readtable(metadataPath,'VariableNamingRule','preserve');

% Create new folder name
fluxFolderName = char(regexp(fbaFolder, '.*[\\/](.*)$', 'tokens', 'once'));

% Create the path to the new folder
newFolder = fullfile(outputFolder, [fluxFolderName '_filtered_samples_in_metadata'] );

% Check if new folder exists
if exist(newFolder,'dir') ~= 7
    mkdir(newFolder)
end

% Prepare fluxes for processing
% Find .mat files in the FBA folder
fbaDirData = what(fbaFolder);

% Rename files for flux processing script
% Remove "microbiota_model_diet_" and ".mat"
resultNames = erase(fbaDirData.mat,{'.mat','FBA_sol_','muWBM_','_male','_female'});

% Get paths to germfree models and remove these paths from the pruned
% database. 
gfIDX = contains(fbaDirData.mat,'gfWBM');

% Remove gf results from paths to prune
resultNames(gfIDX) = [];

% Find results not in the metadata
[~,~,newIndex] = intersect(metadata.ID,resultNames,'stable');

% Do copy the gf samples, so remove gfWBM from ia
resultsToCopy = [fbaDirData.mat(contains(fbaDirData.mat,'gfWBM')); fbaDirData.mat(newIndex+2)];

% Recreate paths for the pruned data
originalFBApaths = append(fbaDirData.path, filesep, resultsToCopy);
prunedPaths = append(newFolder, filesep, resultsToCopy);
originalFBApaths = string(originalFBApaths);
prunedPaths = string(prunedPaths);

% Move files to new folder and rename them
tic
disp('Copy renamed file')
arrayfun(@copyfile, originalFBApaths, prunedPaths)
toc


end