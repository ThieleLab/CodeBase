% run comparison with other resources

mkdir([rootDir filesep 'Comparison_other_GEMs'])

% first extract all CarveMe models form the GitHub repository

% For 90 strains, both a CarveMe and an AGORA2 model exist (see
% matchedModels_CarveMe.xlsx for the list).

modelPath = [rootDir filesep 'Comparison_other_GEMs' filesep 'CarveMe_publishedModels' filesep 'embl_gems'];

% extract models
dInfo = dir(fullfile(modelPath, '**/*.*'));  %get list of files and folders in any subfolder
dInfo = dInfo(~[dInfo.isdir]);
models={dInfo.name};
models=models';
folders={dInfo.folder};
folders=folders';
% remove any files that are not SBML or mat files
delInd=find(~(contains(models(:,1),{'.gz'})));
models(delInd,:)=[];
folders(delInd,:)=[];
for i=1:length(models)
    gunzip([folders{i} filesep models{i}]);
    delete([folders{i} filesep models{i}]);
end

cd([rootDir filesep 'Comparison_other_GEMs']);

% define the resources, models, and info files to test
toTest = {
    'BiGG',[rootDir filesep 'Comparison_other_GEMs' filesep 'BiGG' filesep 'translatedModels'],[rootDir filesep 'Comparison_other_GEMs' filesep 'BiGG' filesep 'mapped_Model_IDs_BiGG.xlsx']
    'CarveMe',[rootDir filesep 'Comparison_other_GEMs' filesep 'CarveMe' filesep 'translatedModels'],[rootDir filesep 'Comparison_other_GEMs' filesep 'CarveMe' filesep 'mapped_Model_IDs_CarveMe.xlsx']
    'gapseq',[rootDir filesep 'Comparison_other_GEMs' filesep 'gapseq' filesep 'translatedModels'],[rootDir filesep 'Comparison_other_GEMs' filesep 'gapseq' filesep 'mapped_Model_IDs_gapseq.xlsx']
    'MAGMA',[rootDir filesep 'Comparison_other_GEMs' filesep 'MAGMA' filesep 'translatedModels'],[rootDir filesep 'Comparison_other_GEMs' filesep 'MAGMA' filesep 'mapped_Model_IDs_MAGMA.xlsx']
    };

for i=1:size(toTest,1)
    cd(toTest{i,1})

    % comparison with Madin,NJC19, and BiGG data
    infoFile = readInputTableForPipeline(toTest{i,3});
    infoFile(:,1:2)=[];

    resourceToComparePath = [rootDir filesep 'Comparison_other_GEMs' filesep toTest{i,1} filesep 'translatedModels'];
    mappedKbasePath = [rootDir filesep 'Comparison_other_GEMs' filesep toTest{i,1} filesep 'kbaseModels'];
    mappedAgora2Path = [rootDir filesep 'Comparison_other_GEMs' filesep toTest{i,1} filesep 'agora2Models'];
    paths = {mappedKbasePath,mappedAgora2Path,resourceToComparePath};
    resources = {'KBase','AGORA2',toTest{i,1}};
    testResourceAgainstExperimentalData(infoFile,paths,resources,numWorkers);
    cd ..
end

%% AGORA2 compared with KBase

cd([rootDir filesep 'Comparison_other_GEMs']);
mkdir('AGORA2')
cd('AGORA2')

infoFile = readInputTableForPipeline('AGORA2_infoFile.xlsx');

paths = {translatedDraftsFolder,refinedFolder};
resources = {'KBase','AGORA2'};
testResourceAgainstExperimentalData(infoFile,paths,resources,numWorkers);
cd ..
