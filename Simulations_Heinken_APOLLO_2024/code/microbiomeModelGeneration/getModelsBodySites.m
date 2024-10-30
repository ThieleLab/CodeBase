
% define path to microbe models
mkdir([rootDir filesep 'Data' filesep 'BodySiteMicrobiomes'])
mkdir(modPath)
modPath=[rootDir filesep 'Data' filesep 'BodySiteMicrobiomes' filesep 'Models'];
mkdir(modPath)

% path to and name of the file with abundance information
abunFilePaths={
    [rootDir filesep 'input' filesep 'normalized_nasal_cavity_abundances.csv']
    [rootDir filesep 'input' filesep 'normalized_vagina_abundances.csv']
    [rootDir filesep 'input' filesep 'normalized_skin_abundances.csv']
    };

for j=1:length(abunFilePaths)
    % get models
    coverage = table2cell(readtable(abunFilePaths{j},'ReadVariableNames',false));
    
    modelsToLoadPath = [rootDir filesep 'data' filesep 'PasolliReconstructions' filesep 'refinedReconstructions'];
    
    for i=2:size(coverage,1)
        if isfile([modelsToLoadPath filesep coverage{i,1} '.mat'])
            copyfile([modelsToLoadPath filesep coverage{i,1} '.mat'],[modPath filesep coverage{i,1} '.mat']);
        end
    end
end

% get the remaining models from AGORA2
modelsToLoadPath = [rootDir filesep 'AGORA2'];
for j=1:length(abunFilePaths)
    % get models
    coverage = table2cell(readtable(abunFilePaths{j},'ReadVariableNames',false));
    for i=2:size(coverage,1)
        if isfile([modelsToLoadPath filesep coverage{i,1} '.mat'])
            copyfile([modelsToLoadPath filesep coverage{i,1} '.mat'],[modPath filesep coverage{i,1} '.mat']);
        end
    end
end
