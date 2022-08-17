
% run comparison with other resources: BiGG, CarveMe, gapseq, MAGMA

cd([rootDir filesep 'ValidationAgainstExperimentalData']);

% define the resources, models, and info files to test
toTest = {'BiGG','CarveMe','gapseq','MAGMA'};

for i=1:length(toTest)
    cd(toTest{i})

    % comparison with Madin,NJC19, and BiGG data
    infoFile = readInputTableForPipeline(['mapped_Model_IDs_' toTest{i} '.xlsx']);
    infoFile(:,1:2)=[];

    resourceToComparePath = [rootDir filesep 'ValidationAgainstExperimentalData' filesep toTest{i} filesep 'translatedModels'];
    paths = {translatedDraftsFolder,refinedFolder,resourceToComparePath};
    resources = {'KBase','AGORA2',toTest{i}};
    testResourceAgainstExperimentalData(infoFile,paths,resources,numWorkers);
    cd ..
end

%% AGORA2 compared with KBase

cd([rootDir filesep 'ValidationAgainstExperimentalData']);
mkdir('AGORA2')
cd('AGORA2')

infoFile = readInputTableForPipeline('AGORA2_infoFile.xlsx');

paths = {translatedDraftsFolder,refinedFolder};
resources = {'KBase','AGORA2'};
testResourceAgainstExperimentalData(infoFile,paths,resources,numWorkers);
cd ..
