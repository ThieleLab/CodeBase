
% path to microbe models
rootDir=pwd;
modPath=[rootDir filesep 'BodySiteMicrobiomes' filesep 'Models'];
mkdir(modPath)

% path to and name of the file with abundance information
abunFilePaths={
   [pwd filesep 'BodySiteMicrobiomes' filesep 'abundances_from_metagenomes' filesep 'normalized_nasal_cavity_abundances.csv']
    [pwd filesep 'BodySiteMicrobiomes' filesep 'abundances_from_metagenomes' filesep 'normalized_vagina_abundances.csv']
    [pwd filesep 'BodySiteMicrobiomes' filesep 'abundances_from_metagenomes' filesep 'normalized_skin_abundances.csv']
    };

for j=1:length(abunFilePaths)
    % get models
    coverage = table2cell(readtable(abunFilePaths{j},'ReadVariableNames',false));
    
    modelsToLoadPath = 'D:\150k_Project\old\curatedReconstructions';
    
    for i=2:size(coverage,1)
        if isfile([modelsToLoadPath filesep coverage{i,1} '.mat'])
            copyfile([modelsToLoadPath filesep coverage{i,1} '.mat'],[modPath filesep coverage{i,1} '.mat']);
        end
    end
end

% get the remaining models from AGORA2
modelsToLoadPath = 'D:\AGORA2_rerunEverything\refinedReconstructions';
for j=1:length(abunFilePaths)
    % get models
    coverage = table2cell(readtable(abunFilePaths{j},'ReadVariableNames',false));
    for i=2:size(coverage,1)
        if isfile([modelsToLoadPath filesep coverage{i,1} '.mat'])
            copyfile([modelsToLoadPath filesep coverage{i,1} '.mat'],[modPath filesep coverage{i,1} '.mat']);
        end
    end
end
