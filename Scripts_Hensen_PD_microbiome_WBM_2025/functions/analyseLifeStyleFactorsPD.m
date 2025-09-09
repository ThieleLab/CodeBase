function [lifeStyleInfluence,filePath] = analyseLifeStyleFactorsPD(paths, confounders, regressions, savePath)
%
% In this function, we will investigate if adding a covariate for any of
% the included lifestyle factors in the metadata will change the regression
% fit. If the covariate is relevant, i.e., p<0.05, for the analysed
% metabolites, the covariate will be added to the analysis. 
%
% INPUTS:
% regressionResults
% preparedInputTable
% preparedMetadata
% rxnsOfInterest
% 
% OUTPUTS:
% newResLifeStyle

lifestyleFactors = {'How_often_do_you_eat_GRAINS',...
    'How_often_do_you_eat_POULTRY_BEEF_PORK_SEAFOOD_EGGS',...
    'How_often_do_you_eat_FRUITS_or_VEGETABLES',...
    'How_often_do_you_eat_NUTS',...
    'How_often_do_you_eat_YOGURT',...
    'Do_you_smoke',...
    'Do_you_drink_caffeinated_beverages'};

% Next, we will add each lifestyle variable to the list of confounders and
% investigate the influence of this confounder on the regression fit. 
newResLifeStyle = cell(length(lifestyleFactors),1); % Preallocate cell array
for i=1:length(lifestyleFactors)
    i
    confoundersLF = [confounders,lifestyleFactors(i)]; % Append covariate
    %try
    [~,~,~, ~, altRegressions] = performParkinsonAnalysis(paths.fluxPath,paths.metadataPath, confoundersLF, paths.parkinsonFluxes);

    % Perform a likelihood-ratio test by performing the chi-square for the
    % difference in model deviance from the nested model (mdl1) to the expanded
    % model (mdl2).
    lrt = @(mdl1,mdl2) chi2cdf(mdl1.Deviance-mdl2.Deviance, mdl1.DFE-mdl2.DFE,'upper');
    % Append model deviances and the change in deviance
    devLrt = @(mdl1,mdl2) [mdl1.Deviance, mdl2.Deviance, mdl1.Deviance-mdl2.Deviance, mdl1.DFE-mdl2.DFE, ...
        mdl1.ModelCriterion.AIC-mdl2.ModelCriterion.AIC,mdl1.ModelCriterion.BIC-mdl2.ModelCriterion.BIC...
        lrt(mdl1,mdl2)]; 
    
    % Iterate over each field in the regression structures
    fNames = string(fieldnames(regressions));
    modDiff = arrayfun(@(x) devLrt(regressions.(x), altRegressions.(x)), fNames,'UniformOutput',false);
    modDiff = array2table(vertcat(modDiff{:}),'VariableNames',...
        {'Deviance_1','Deviance_2','delta_Deviance','delta_DF','delta_AIC','delta_BIC','chisq_p'});
    modDiff = addvars(modDiff,fNames,'NewVariableNames','Reaction','Before',1);
    modDiff = addvars(modDiff,repmat(string(lifestyleFactors(i)),height(modDiff),1),'NewVariableNames','Variable','After','Reaction');
    modDiff = sortrows(modDiff,'chisq_p','ascend');
    newResLifeStyle{i} = modDiff;
    %end
end

lifeStyleInfluence = vertcat(newResLifeStyle{:});
lifeStyleInfluence = sortrows(lifeStyleInfluence,'chisq_p','ascend');
lifeStyleInfluence.Reaction = replace(lifeStyleInfluence.Reaction,'_bc_','[bc]');

lifeStyleInfluence = lifeStyleInfluence(matches(lifeStyleInfluence.Reaction,paths.rxnsOfInterest),:);

lifeStyleInfluence.Reaction = renameVmhToMetName(lifeStyleInfluence.Reaction);

% G = groupsummary(lifeStyleInfluence,"Variable","mean","chisq_p")

% Save table
filePath = fullfile(savePath,'lifeStyleCovariateAnalysis.csv');
writetable(lifeStyleInfluence,filePath)


%%

% Next, we will add each lifestyle variable to the list of confounders and
% investigate the influence of this confounder on the regression fit. 
% newResLifeStyle = cell(length(lifestyleFactors),1); % Preallocate cell array
% for i=1:length(lifestyleFactors)
% 
%     confoundersLF = [confounders,lif
% end
% 
% estyleFactors(i)]; % Append covariate
%     regFormula = strcat('Case_status~','Flux+',strjoin(confoundersLF,'+')); % Define regression formula
%     newResults = performRegressions(preparedInputTableNew,preparedMetadataNew,regFormula); % Perform regressions
% 
%     % Filter on regression coefficients of the added covariates and
%     % generate a single table with all results
%     newResults = rmfield(newResults,fieldnames(regressionResults)); 
%     newResults = struct2cell(newResults);
%     newResults = vertcat(newResults{:});
%     newResults = sortrows(newResults,'pValue','ascend');
%     newResults.Formula = string(newResults.Formula);
%     newResLifeStyle{i} = newResults;
% end
% 
% % Concatenate all regression results in a single table
% newResLifeStyle = vertcat(newResLifeStyle{:});
% 
% newResLifeStyle.Reaction = renameVmhToMetName(newResLifeStyle.Reaction);
% 
% % Conclusion: After FDR correction, none of the lifestyle factors
% % influenced the regression model. Thus, the no changes will be made on the
% % analysed list of confounders. 
% 
% % Save table
% filePath = fullfile(savePath,'lifeStyleCovariateAnalysis.csv');
% writetable(newResLifeStyle,filePath)