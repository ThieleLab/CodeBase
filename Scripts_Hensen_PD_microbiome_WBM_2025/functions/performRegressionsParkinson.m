function [results,regressions] = performRegressionsADRC(data,metadata,formula,regressionType,responseOrdination, exponentiateLogOdds)
% This functions performs linear or logistic regressions on either flux
% data or microbial relative abundances. The function supports control
% variables & moderators in the regression formulas. 
%
% USAGE:
%    results = performRegressions(data,metadata,formula)
% 
% INPUT
% data:     Table with flux or microbiome data. The first column must be the ID
%           column.
% metadata: Metadata table. The first column must be the ID column.
% formula: Must contain either "Flux or "relative_abundance" as predictor.
% exponentiateLogOdds: Logical variable on if a logistic regression estimate
% should be exponentiated to obtain the odds ratios.
%
% .. Author:
%       - Tim Hensen, 2024

%regressionType = 'multinomial';

if strcmp(regressionType,'ordinal_multinomial') & nargin<5
    error('Please give response variable ordination')
end

if nargin < 6
    exponentiateLogOdds = false;
end

% Check if formula is correct
if contains(formula,'Flux')
    value = 'Flux';
    name = 'Reaction';
elseif contains(formula,'relative_abundance')
    value = 'relative_abundance';
    name = 'Taxa';
else
    error('Please add correct formula input')
end

% Get names of metabolites or taxa
names = data.Properties.VariableNames;
names(matches(names,'ID'))=[];

% Stack
data = stack(data,names,'NewDataVariableName',value,'IndexVariableName',name);
data.ID = strrep(data.ID,'HM_','');
% combine with metadata
data = innerjoin(data,metadata,'Keys','ID');

% Filter on variables included in the formula

% Find response and predictors 
includedVars = strtrim(strsplit(formula,'[+\|*]','DelimiterType','RegularExpression'));
variables = [strtrim(strsplit(includedVars{1},'~')) includedVars(2:end)];
response = variables(1);
predictors = variables(2:end);
% Remove moderator term (defined by use of ":") from predictors
predictors = predictors(~contains(predictors, ":"));

% Check if regression should be linear or logistic

switch regressionType
    case 'linear'
        responseDistribution = 'normal';
    case 'logistic'
        responseDistribution = 'binomial';
        if ~isempty(responseOrdination)
            data.(string(response)) = grp2idx( categorical(data.(string(response)),responseOrdination) )-1;
        else
            data.(string(response)) = grp2idx( categorical(data.(string(response))) )-1;
        end
            
    case 'multinomial'
        data.(string(response)) = categorical(data.(string(response)));
        responseDistribution = 'multinomial';
    case 'ordinal_multinomial'
        data.(string(response)) = categorical(data.(string(response)), responseOrdination);
        responseDistribution = 'multinomial';
end

respVarData = data.(string(response));
if isscalar(unique(respVarData))
    error('Make sure that that there is more than one category in the response variable.')
end

% Filter on variables of interest
includedVars = [{name} response predictors];
data = data(:,includedVars);

% Find metabolite groups
[groups,groupnames] = findgroups(data.(name));
groupnames = string(groupnames);
fieldNames = matlab.lang.makeValidName(groupnames);

% Check current matlab version. Perform firth regression for improved
% robustness if the matlab version is above 2024b.

% Find current matlab version
matlabVersion = ver;
matlabVersion = str2double(string({matlabVersion.Version}'));

firthRegression = false;
if matlabVersion(1)>=24.2 && ~contains(regressionType,'multinomial') % Check if the current matlab version is greater than 2024a
    firthRegression = true;
end

% Perform regressions
regressions = struct;
for i = 1:length(unique(groups))
    lastwarn('')

    if firthRegression == true
        mdl = fitglm(data(groups==i,:),formula,'Distribution',responseDistribution,'LikelihoodPenalty','jeffreys-prior'); % Firth's regression
    elseif strcmp(regressionType,'ordinal_multinomial')
        mdl = fitmnr(data(groups==i,:),formula,'ModelType','ordinal','Link','probit');
    elseif strcmp(regressionType,'multinomial')
        mdl = fitmnr(data(groups==i,:),formula,'Link','probit');
    else
        mdl = fitglm(data(groups==i,:),formula,'Distribution',responseDistribution);
    end     

    % If no good fit could be found, do not save the result
    warnMsg = lastwarn;
    if ~isempty(warnMsg)
        regressions.(fieldNames(i)) = {};
    else
        regressions.(fieldNames(i)) = mdl;
    end
end

% Assuming a logistic regression, create the following table

% Create empty table
varNames = {name,'Formula','Predictor','Regression type','N','estimate','low','high','SE','tStat','pValue','FDR','R2'};
varTypes = [repmat({'string'},1,4),repmat({'double'},1,length(varNames)-4)];
generalTable = table('Size',[length(groupnames),length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);

% Prefill table
% formula4table = [mdl.Formula.ResponseName, '~', mdl.Formula.LinearPredictor];
generalTable.Formula = repmat(string(formula),length(generalTable.(string(name))),1);
generalTable.('Regression type') = repmat(string(responseDistribution),length(generalTable.(string(name))),1);

% Find the first non empty field in the regressions variable
nonEmptyResults = find(~structfun(@isempty,regressions));

% Preallocate structure
results = struct;
predictor = 'NotDefined';

% Get the number of predictors in the data
if ~isempty(nonEmptyResults)
    predictors = regressions.(fieldNames(nonEmptyResults(1))).CoefficientNames(2:end);


    for j=2:length(predictors)+1
        % Populate tables 
        predictorTable = generalTable;
        for i=1:length(fieldNames)
            mdl = regressions.(fieldNames(i));
            % Check if a result could be found
            if ~isempty(mdl)
                % Find predictor name
                predictor = string(matlab.lang.makeValidName(mdl.CoefficientNames(j)));
            
                % Add feature name
                predictorTable.(name)(i) = string(mdl.Variables.(name)(1));
                
                % Add sample number
                predictorTable.N(i) = mdl.NumObservations;
                % Add regression estimate
                predictorTable{i,{'estimate','SE','tStat','pValue'}} = mdl.Coefficients{j,:};
                % Add 95% CI for predictor
                ci=coefCI(mdl);
                predictorTable{i,{'low','high'}} = ci(j,:);
                % Calculate odds ratios if responseDistribution = 'binomial'
                if matches(responseDistribution,'binomial')
                    if exponentiateLogOdds
                        predictorTable{i,{'estimate','low','high'}} = exp(predictorTable{i,{'estimate','low','high'}});
                    else
                        predictorTable{i,{'estimate','low','high'}} = predictorTable{i,{'estimate','low','high'}};
                    end
                end
    
                % Add adjusted R2
                predictorTable.R2(i) = mdl.Rsquared.Adjusted;
            else
                % Make the missing data nans
                predictorTable{i,{'N','estimate','low','high','SE','tStat','pValue','FDR','R2'}}=nan;
                predictorTable.(name)(i) = fieldNames(i);
            end
        end
        
        % Add FDR values 
        % predictorTable.FDR = mafdr(predictorTable.pValue,'BHFDR',true);
        predictorTable.FDR = fdrBHadjustment(predictorTable.pValue); % Local alternative. 

        % Sort by significance
        predictorTable = sortrows(predictorTable,'pValue','ascend');
    
        % Add predictor
        predictorTable.Predictor = repmat(predictor,length(predictorTable.(string(name))),1);
    
        % Add result to structure
        results.(predictor) = predictorTable;
    end

else
    results.(predictor) = {'No results'};
    warning('None of the variables could be investigated.')
end
end