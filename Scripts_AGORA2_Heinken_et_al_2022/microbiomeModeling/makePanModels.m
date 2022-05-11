
% AGORA2
allpanFolder=[rootDir filesep 'PanModelsAGORA2'];
mkdir(allpanFolder)
cd(allpanFolder)

taxa={
    'Species'
    'Genus'
    'Family'
    'Order'
    'Class'
    'Phylum'
    };

for i=1:length(taxa)
    panFolder=[allpanFolder filesep taxa{i}];
    mkdir(panFolder)
    createPanModels(refinedFolder,panFolder,taxa{i},numWorkers,'AGORA2_infoFile.xlsx');
    testResultsFolder=[allpanFolder filesep 'TestResults_' taxa{i}];
    [notGrowing,biomassFluxes] = plotBiomassTestResults(panFolder, taxa{i}, 'numWorkers',numWorkers, 'testResultsFolder', testResultsFolder);
    [tooHighATP,atpFluxes] = plotATPTestResults(panFolder, taxa{i}, 'numWorkers',numWorkers, 'testResultsFolder', testResultsFolder);
end

cd(rootDir)

