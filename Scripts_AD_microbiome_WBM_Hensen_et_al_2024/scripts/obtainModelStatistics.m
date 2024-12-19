function modelStatPath = obtainModelStatistics(workingDir, fluxResultPath, hmDir, microbiomeDir)
% This function will load the HM and microbiome models in the metadata and
% extract for the HM models: reaction count, metabolite count, constraint
% count, for male and female samples. For the microbiome models, the number
% of reactions and number of subsystems will be extracted. Then, we will
% calculate the mean and standard deviations for these metrics and save the
% output. 
%
% USAGE:
%    modelStatPath = obtainModelStatistics(workingDir, fluxResultPath, hmDir, microbiomeDir)
%
% INPUTS:
%    workingDir:        (char) Path to directory with WBM and microbiome models
%    fluxResultPath:    (char) Path to file with the processed flux results
%                       of 1065 samples
%    hmDir:             (char) Path to directory with WBM models
%    microbiomeDir:     (char) Path to directory with microbiome community models
%
% OUTPUTS:
%    modelStatPath: (char) Path to excel file with results
%
% NOTE:
%    This function should be run after fully running AnalysisPipeline.Rmd
%
% .. Author: Tim Hensen (12/2024)

% Set path to model statistics file
modelStatPath = [workingDir filesep 'modelStats.xlsx'];

% Get model names
processNames = @(x) string(erase(x,{'.mat','microbiota_model_diet_Sample'}));
hmNames = processNames(what(hmDir).mat);
microbiomeNames = processNames(what(microbiomeDir).mat);

% Get paths to models
hmPaths = string(append(what(hmDir).path, filesep, what(hmDir).mat));
microbiomePaths = string(append(what(microbiomeDir).path, filesep, what(microbiomeDir).mat));


% Load metadata IDs 
fluxes = readtable(fluxResultPath,'VariableNamingRule','preserve');
IDs = string(fluxes.ID);

% Remove models not in the metadata
[~,hmPathsToRemove] = setdiff(hmNames,IDs);
hmPaths(hmPathsToRemove)=[];
[~,microbiomePathsToRemove] = setdiff(microbiomeNames,IDs);
microbiomePaths(microbiomePathsToRemove)=[];

% Get ID and sex information
IDs = fluxes.ID;
Sex = fluxes.Sex;

% Preallocate variables for model statistics
wrxns = nan(length(hmPaths),1);
wmets = nan(length(hmPaths),1);
wconstraints = nan(length(hmPaths),1);
mrxns = nan(length(hmPaths),1);
msubsystems = nan(length(hmPaths),1);

% Set parallel pool
numWorkers = feature('numCores');
if numWorkers > 1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end

% Start for loop
parfor i=1:length(hmPaths)
    disp(i)
    % Load models
    Wbm = load(hmPaths(i)); 
    Wbm = Wbm.(string(fieldnames(Wbm)));
    microbiome = load(microbiomePaths(i)); 
    microbiome = microbiome.(string(fieldnames(microbiome)));

    % Calculate the necessary statistics
    wrxns(i) = numel(Wbm.rxns);
    wmets(i) = numel(Wbm.mets);
    wconstraints(i) = numel(Wbm.d);
    mrxns(i) = numel(microbiome.rxns);
    subsystems = unique(string(microbiome.subSystems));
    subsystems = subsystems(~matches(subsystems,""));
    msubsystems(i) = numel(subsystems);
end

% Organise data to table
modelStats = table(IDs, Sex, wrxns, wmets, wconstraints, mrxns, msubsystems);
modelStats.Properties.VariableNames = ...
    {'ID','sex','WBM_reactions','WBM_metabolites','WBM_constraints',...
    'Microbiome_reactions','Microbiome_subsystems'};

% Calculate summary statistics per sex

femaleRows = matches(modelStats.sex,'Female');

% Obtain summary statistics
summaryStats = zeros(5,4);
% Calculate mean values for female samples 
summaryStats(:,1) = table2array(varfun(@mean,modelStats(femaleRows,3:end)))';
% Calculate SD values for female samples 
summaryStats(:,2) = table2array(varfun(@std,modelStats(femaleRows,3:end)))';
% Calculate mean values for male samples 
summaryStats(:,3) = table2array(varfun(@mean,modelStats(~femaleRows,3:end)))';
% Calculate SD values for male samples 
summaryStats(:,4) = table2array(varfun(@std,modelStats(~femaleRows,3:end)))';

% Convert to table
summaryStats = array2table(summaryStats,...
    "RowNames",modelStats.Properties.VariableNames(3:end)',...
    "VariableNames",{'Mean female samples','SD female samples','Mean male samples','SD male samples'});

% Save table
writetable(modelStats, modelStatPath,'Sheet','Data')

% Add sheet to table with WBM statistics
writetable(summaryStats,modelStatPath,'Sheet','Statistics','WriteRowNames',true);
end