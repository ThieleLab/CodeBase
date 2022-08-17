
% Cluster reaction presence in subsets of AGORA2 models 

currentDir=pwd;
cd([rootDir filesep 'modelProperties'])

% define certain classes for testing
toTest={
    'Bacilli'
    'Gammaproteobacteria'
    };

for i=1:length(toTest)
    % refined reconstructions
    mkdir([toTest{i} '_refined'])
    cd([toTest{i} '_refined'])
    % get subset
    [extractedSubset,subsetFolder] = extractReconstructionResourceSubset(refinedFolder, infoFilePath, 'Class', toTest{i}, [pwd filesep 'Reconstructions']);
    % determine reaction presence in subset
    
    getReactionMetabolitePresence(subsetFolder,pwd,reconVersion,numWorkers)
    % cluster subset by reaction presence
    producetSNEPlots(pwd,infoFilePath,reconVersion)
    
    rmdir(subsetFolder,'s')
    cd ..
    
     % draft reconstructions
    mkdir([toTest{i} '_draft'])
    cd([toTest{i} '_draft'])
    % get subset
    [extractedSubset,subsetFolder] = extractReconstructionResourceSubset(translatedDraftsFolder, infoFilePath, 'Class', toTest{i}, [pwd filesep 'Reconstructions']);
    % determine reaction presence in subset
    
    getReactionMetabolitePresence(subsetFolder,pwd,reconVersion,numWorkers)
    % cluster subset by reaction presence
    producetSNEPlots(pwd,infoFilePath,reconVersion)
    
    rmdir(subsetFolder,'s')
    cd ..
end

cd(currentDir)
