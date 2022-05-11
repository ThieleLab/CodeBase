
% Cluster reaction presence in subsets of AGORA2 models 

currentDir=pwd;
cd(propertiesFolder)

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
    
    getReactionMetabolitePresence(subsetFolder,[pwd filesep 'ReactionPresence'],reconVersion,numWorkers)
    % cluster subset by reaction presence
    producetSNEPlots([pwd filesep 'ReactionPresence'],infoFilePath,reconVersion)
    
    rmdir(subsetFolder,'s')
    
     % draft reconstructions
    mkdir([toTest{i} '_draft'])
    cd([toTest{i} '_draft'])
    % get subset
    [extractedSubset,subsetFolder] = extractReconstructionResourceSubset(translatedDraftsFolder, infoFilePath, 'Class', toTest{i}, [pwd filesep 'Reconstructions']);
    % determine reaction presence in subset
    
    getReactionMetabolitePresence(subsetFolder,[pwd filesep 'ReactionPresence'],reconVersion,numWorkers)
    % cluster subset by reaction presence
    producetSNEPlots([pwd filesep 'ReactionPresence'],infoFilePath,reconVersion)
    
    rmdir(subsetFolder,'s')
    cd ..
end

cd(currentDir)
