
clear all
rootDir = pwd;

% define path to microbe models
modPath=[rootDir filesep 'Data' filesep 'GutMicrobiomes' filesep 'Models'];
mkdir(modPath)

% get models
coverage = table2cell(readtable([rootDir filesep 'data' filesep 'GutMicrobiomes' filesep 'normalizedCoverage.csv'],'ReadVariableNames',false));

modelsToLoadPath = [rootDir filesep 'data' filesep 'AlmeidaReconstructions' filesep 'refinedReconstructions'];

for i=2:size(coverage,1)
    if isfile([modelsToLoadPath filesep coverage{i,1} '.mat'])
        copyfile([modelsToLoadPath filesep coverage{i,1} '.mat'],[modPath filesep coverage{i,1} '.mat']);
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
