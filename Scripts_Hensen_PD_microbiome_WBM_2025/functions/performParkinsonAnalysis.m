function [results,preparedInputTable,preparedMetadata, preparedRawInputTable] = performParkinsonAnalysis(dataPath, metadataPath, confounders, savePath)
% This function prepares the flux/reads tables and then performs
% regressions.
%
% AUTHOR: Tim Hensen, January 2025


% Load input data and metadata
inputTable = readtable(dataPath,'VariableNamingRule','preserve');
inputTable.ID = erase(inputTable.ID,{'muWBM_','_male','_female'});
metadata = readtable(metadataPath,'VariableNamingRule','preserve');


%dataForTim = outerjoin(metadata,inputTable,'MergeKeys',true,'Type','left');
% Save results
%writetable(dataForTim,fullfile(saveDir,'processedPDdataForTim.csv'));


% Remove sex information if present
if any(matches(inputTable.Properties.VariableNames,'sex','IgnoreCase',true))
    inputTable.Sex = [];
end

% Make sure that the inputTable and metadata have the same samples
[~,idxa,idxb] = intersect(string(metadata.ID),string(inputTable.ID),'stable');

if numel(idxb)<length(metadata.ID)
    error('COBRA:BadInput', 'No overlapping samples could be found between the reads/flux table and the metadata table.')
else
    metadata = metadata(idxa,:);
    inputTable = inputTable(idxb,:);
end

preparedRawInputTable = inputTable;

% Before performing a log transformation, check if any of the values are
% negative. If so, 1) set all zeros to nan, and 2) add the minimum value
% that makes all values positive
inputMatrix = table2array(inputTable(:,2:end));
if any(any(inputMatrix<0))
    inputMatrix(inputMatrix==0) = nan;
    inputMatrix = inputMatrix + (abs(min(min(inputMatrix))) *1.1);
end

% Transform the inputData by performing a log2 transformation and
% z-scaling.
log2Data = log2(inputMatrix);
log2Data(isinf(log2Data)) = nan;
inputTable{:,2:end}=normalize(log2Data);
inputTable.ID = erase(inputTable.ID,{'muWBM_','_female','_male'});
preparedInputTable = inputTable;

% Transform the metadata, age, bmi, and total  sequence count for 
% statistical analysis.
metadata.Age_at_collection = normalize(metadata.Age_at_collection);
metadata.BMI = normalize(log10(metadata.BMI));
metadata.total_sequences = normalize(log10(metadata.total_sequences));
preparedMetadata = metadata;

% % Set confounding variables
% confounders = {'Age_at_collection','Sex','Age_at_collection',...
%     'Do_you_drink_alcohol','Laxatives','Probiotic','Antihistamines',...
%     'Depression_anxiety_mood_med','Pain_med','Sleep_aid','total_sequences','Antibiotics_current'};

% Prepare regression formula
regFormula = strcat('Case_status~','Flux+',strjoin(confounders,'+'));

% Perform regressions
[results,regressions] = performRegressions(preparedInputTable,preparedMetadata,regFormula);

% Add PD and control samples to the flux results
results = addPdSampleCounts(results,regressions);

% If present, remove results from carnitine and O-Octadecanoyl-R-Carnitine
results.Flux(matches(results.Flux.Reaction,{'DM_stcrn[bc]','DM_crn[bc]'}),:)=[];

% Add the number of PD samples and controls to the regression results

% Save results
writetable(results.Flux,savePath)

end


function results = addPdSampleCounts(results,regressions)
% Find the number of PD samples and controls in the regressions

% Prepare table with reactions and sample numbers
testedRxns = results.Flux.Reaction;
sampNumbers = table(testedRxns,zeros(numel(testedRxns),1),zeros(numel(testedRxns),1),'VariableNames',{'Reaction','PD_N','Control_N'});

% Get the variable info for the regression data of reaction i
for i=1:numel(testedRxns)
    rxnToFill = matlab.lang.makeValidName(sampNumbers.Reaction(i));

    % Get data
    data = regressions.(rxnToFill).Variables;

    % Remove nan values in the flux data
    data(isnan(data.Flux),:) = [];

    % Add the number of cases and controls to the sampNumbers table
    sampNumbers.PD_N(i) = sum(data.Case_status==1);
    sampNumbers.Control_N(i) = sum(data.Case_status==0);
end

% Add data to results
results.Flux  = outerjoin(results.Flux,sampNumbers,'Keys','Reaction','Type','left','MergeKeys',true);
results.Flux = movevars(results.Flux,'PD_N','After','N');
results.Flux = movevars(results.Flux,'Control_N','After','PD_N');

% Sort results again
results.Flux = sortrows(results.Flux,'pValue','ascend');
end