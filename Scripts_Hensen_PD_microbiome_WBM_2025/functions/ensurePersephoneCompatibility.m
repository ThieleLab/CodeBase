function OK = ensurePersephoneCompatibility(fluxPath)
% Rename the solution field speciesBIO to taxonNames in the flux .mat files

% Get flux solutions
fbaDirData = what(fluxPath).mat;

% Remove gf results from paths to prune
gfIDX = contains(fbaDirData,'gfWBM');
fbaDirData(gfIDX) = [];

% Get paths to FBA solutions
solPaths = string(fullfile(fluxPath,fbaDirData));

% Define field names to change
for i=1:length(solPaths)
    taxonNames = load(solPaths(i)).speciesBIO; % Load old name
    save(solPaths(i),'taxonNames','-append'); % Append new name to results
end

OK = true;
end