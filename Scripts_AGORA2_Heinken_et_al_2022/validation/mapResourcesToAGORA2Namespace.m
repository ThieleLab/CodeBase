
% map BiGG, CarveMe, gapseq, and MAGMA to VMH namespace to enable
% comparison with AGORA2 against independent experimComparison_other_GEMsental data

% define the resources, folders, and info files to test
toTest = {
    'BiGG',[rootDir filesep 'Comparison_other_GEMs' filesep 'BiGG' filesep 'models'],[rootDir filesep 'Comparison_other_GEMs' filesep 'BiGG' filesep 'mapped_Model_IDs_BiGG.xlsx']
    'CarveMe',[rootDir filesep 'Comparison_other_GEMs' filesep 'CarveMe' filesep 'embl_gems'],[rootDir filesep 'Comparison_other_GEMs' filesep 'CarveMe' filesep 'mapped_Model_IDs_CarveMe.xlsx']
    'gapseq',[rootDir filesep 'Comparison_other_GEMs' filesep 'gapseq' filesep 'models'],[rootDir filesep 'Comparison_other_GEMs' filesep 'gapseq' filesep 'mapped_Model_IDs_gapseq.xlsx']
    'MAGMA',[rootDir filesep 'Comparison_other_GEMs' filesep 'MAGMA' filesep 'models'],[rootDir filesep 'Comparison_other_GEMs' filesep 'MAGMA' filesep 'mapped_Model_IDs_MAGMA.xlsx']
    };

cd('Comparison_other_GEMs')
for i=1:size(toTest,1)
    cd(toTest{i,1})
     mapResourceToAGORA2(toTest{i,2},translatedDraftsFolder,refinedFolder,toTest{i,3},toTest{i,1})
     cd ..
end
cd ..
