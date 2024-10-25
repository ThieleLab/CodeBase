
% initialize the COBRA Toolbox and solvers
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

% path to microbe models
modPath=[pwd filesep 'Models'];

numWorkers=20;

dietFilePath=[pwd filesep 'AverageEuropeanDiet_modified'];

computeProfiles=false;

%% split into multiple runs for better performance

% get split coverages files
dInfo = dir('split_coverages');
fileList={dInfo.name};
fileList=fileList';
fileList(find(strncmp(fileList,'.',1)),:)=[];
[~, reindex] = sort(str2double(regexp(fileList, '\d+', 'match', 'once' )));
fileList = fileList(reindex);

for i=1:length(fileList)
    % path to the abundance file
    abunFilePath=[pwd filesep 'split_coverages' filesep fileList{i}];
    
    % path where to save results
    resPath=[pwd filesep 'MicrobiomeModels' filesep 'MicrobiomeModels_' num2str(i)];
    
    [init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] = initMgPipe(modPath, abunFilePath, computeProfiles, 'dietFilePath', dietFilePath, 'resPath', resPath, 'numWorkers', numWorkers);   
end
