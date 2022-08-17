
% map BiGG, CarveMe, gapseq, and MAGMA to VMH namespace to enable
% comparison with AGORA2 against independent experimental data

% define the resources, folders, and info files to test
toTest = {
    'BiGG',[rootDir filesep 'ValidationAgainstExperimentalData' filesep 'BiGG' filesep 'models'],[rootDir filesep 'ValidationAgainstExperimentalData' filesep 'BiGG' filesep 'mapped_Model_IDs_BiGG.xlsx']
    'CarveMe',[rootDir filesep 'ValidationAgainstExperimentalData' filesep 'CarveMe' filesep 'models'],[rootDir filesep 'ValidationAgainstExperimentalData' filesep 'CarveMe' filesep 'mapped_Model_IDs_CarveMe.xlsx']
    'gapseq',[rootDir filesep 'ValidationAgainstExperimentalData' filesep 'gapseq' filesep 'models'],[rootDir filesep 'ValidationAgainstExperimentalData' filesep 'gapseq' filesep 'mapped_Model_IDs_gapseq.xlsx']
    'MAGMA',[rootDir filesep 'ValidationAgainstExperimentalData' filesep 'MAGMA' filesep 'models'],[rootDir filesep 'ValidationAgainstExperimentalData' filesep 'MAGMA' filesep 'mapped_Model_IDs_MAGMA.xlsx']
    };

cd('ValidationAgainstExperimentalData')
for i=1:size(toTest,1)
    cd(toTest{i,1})
     mapResourceToAGORA2(toTest{i,2},translatedDraftsFolder,refinedFolder,toTest{i,3},toTest{i,1})
     cd ..
end
cd ..
