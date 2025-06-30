function [fullModRes, filePath] = moderationAnalysisPD(preparedInputTable, preparedMetadata, confounders, rxnsOfInterest, saveDir)
% Investigate the influence of age and sex on the associations between the
% predicted fluxes and PD status. 
%
% INPUT:
% preparedInputTable
% preparedMetadata
% saveDir
%
% OUTPUT:
% fullModRes
% filePath

% Filter flux table on reactions of interest
preparedInputTable = preparedInputTable(:,["ID",rxnsOfInterest']);

% Define parameters for moderation analysis
regressionParam.response = 'Case_status'; % Response variable
regressionParam.confounders = confounders; % Control variables
regressionParam.interactionPvalThreshold = 1; % Include all results

% Define age as the interaction variable of interest flux:age
regressionParam.modVar = 'Age_at_collection';

% Perform interaction regressions on age 
ageModRes = runModerationAnalysis(preparedInputTable, preparedMetadata,regressionParam);
%%
% Set Sex as the next interaction variable
regressionParam.modVar = 'Sex'; 

% Perform moderation analysis on Sex 
sexModRes = runModerationAnalysis(preparedInputTable, preparedMetadata,regressionParam);
%%
% Filter on full
sexModRes = sexModRes(matches(sexModRes.Cohort,'Full'),:);
sexModRes = removevars(sexModRes,'Cohort');

%%

% Prepare metadtata table for investigation on the influence of dietary
% habits on the flux-PD associations.
% metadata = preparedMetadata;
% 
% ffq = {'How_often_do_you_eat_GRAINS',...
%     'How_often_do_you_eat_POULTRY_BEEF_PORK_SEAFOOD_EGGS',...
%     'How_often_do_you_eat_FRUITS_or_VEGETABLES',...
%     'How_often_do_you_eat_NUTS',...
%     'How_often_do_you_eat_YOGURT'};
% baseCat = {'At least once a day','At least once a day','At least once a day','At least once a day','Few times a week'};
% regressionParam.interactionPvalThreshold = 0.05; % Include all results
% 
% modResDiet = cell(length(ffq),1);
% for i=1:length(ffq)
%     metadata.(ffq{i}) = string(matches(metadata.(ffq{i}), baseCat{i}));
%     % Perform moderation analysis on Sex 
%     regressionParam.modVar = ffq{i}; 
%     modResDiet{i} = runModerationAnalysis(preparedInputTable, metadata,regressionParam);
% end

% Set Sex as the next interaction variable



%%

% Combine interaction results
fullModRes = [ageModRes;sexModRes];

% Save results
filePath = fullfile(saveDir,'fluxAgeSexModerationStats.csv');
writetable(fullModRes,filePath)

end



function stratResTable = runModerationAnalysis(preparedInputTable, preparedMetadata,regressionParam)


% Define parameters for moderation analysis
response = regressionParam.response;
confounders = regressionParam.confounders;
modVar = regressionParam.modVar;
interactionPvalThreshold = regressionParam.interactionPvalThreshold;

% If the modVar is already in the confounder list, remove this variable
% from the confounder list.
confounders(matches(confounders,modVar))=[];


% Define regression formula
regFormula = strcat(response,'~','Flux+',modVar,'+',strjoin(confounders,'+'),'+Flux:',modVar);

% Make sure that the moderation variable is not numeric
proceed = true;
if isnumeric(preparedMetadata.(modVar)) & numel(unique(preparedMetadata.(modVar)))>2
    proceed = false;
    % modvarDat = cellstr(num2str(preparedMetadata.(modVar)));
    % modvarDat(matches(modvarDat,{'0'})) = {'No'};
    % modvarDat(matches(modvarDat,{'1'})) = {'Yes'};
    % preparedMetadata.(modVar) = modvarDat;
end

% Perform regression with interaction effect

results = performRegressions(preparedInputTable,preparedMetadata,regFormula);

% Extract interaction effects 
resFields = string(fieldnames(results));
intResults = results.(resFields(end));
intResultsOri = intResults;

% Test if any further interaction analyses will be performed

if ~any(matches(resFields,'NotDefined')) & proceed == true
    proceed = true;
end

if proceed == true
    % Filter on reactions of interest
    intResults( isnan(intResults.estimate) | intResults.pValue > interactionPvalThreshold,:) = [];
    rxnsToTest = intResults.Reaction';
    % Test if any further interaction analyses will be performed
    proceed = false;
    if ~isempty(rxnsToTest)
        proceed = true;
    end
end

if proceed == true
    % Prune flux table
    interactionInputTable = preparedInputTable(:,["ID" rxnsToTest]);
    
    % Find binary groups in interaction variable
    [~,cohort] = findgroups(preparedMetadata.(modVar));

    
    % Stratify metadata based on group status
    metadataStrat = cellfun(@(x) preparedMetadata(matches(preparedMetadata.(modVar), x), :), cohort, 'UniformOutput', false);
    
    % Redefine regression formula
    regFormula = strcat(response,'~','Flux+',strjoin(confounders,'+')); 

    % Perform regressions on both strata
    stratRes = cellfun(@(md) performRegressions(interactionInputTable, md, regFormula),metadataStrat, 'UniformOutput', false);

    % Remove cell entry if no regressions could be performed on a stratum
    rmRes = cellfun(@(x) any(matches(fieldnames(x),{'NotDefined'})), stratRes,'UniformOutput',true);
    stratRes = stratRes(~rmRes);
    cohort = cohort(~rmRes);

    if isempty(stratRes)
        proceed = false;
    end
end

if proceed == true

    % Unnest stratified results and add cohort label
    addCohortFun = @(x,y) addvars(x.Flux,repmat(string(y),height(x.Flux),1), 'Before','N','NewVariableNames','Cohort');
    stratRes = cellfun(addCohortFun,stratRes,cohort,'UniformOutput',false);

    stratRes = cellfun(@(x) convertvars(x,'Formula','string'),stratRes,'UniformOutput',false);
    
    % Also add cohort label to the interaction results
    stratRes{3} = addvars(intResults,repmat("Full",height(intResults),1), 'Before','N','NewVariableNames','Cohort');
    
    % Create a single table for the moderation analysis results
    stratResTable = vertcat(stratRes{:});
    stratResTable = sortrows(stratResTable,'Reaction');
else
    stratResTable = intResultsOri;
    disp('No interaction effect found for any of the reactions and moderation variables')
end

end

