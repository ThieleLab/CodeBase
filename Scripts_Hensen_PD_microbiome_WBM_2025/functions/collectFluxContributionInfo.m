function [mContributions,filteredfmLinksSummary] = collectFluxContributionInfo(mContributionTablePath,filteredfmLinksSummaryPath)

% Load flux metabolite contribution table
mContributions = readtable(mContributionTablePath,'VariableNamingRule','preserve');

% Rename variables for supplementary tables
currNames = mContributions.Properties.VariableNames;
newNames = {'WBM metabolite in blood','Associated microbial species','Bootstrapped mean','2.5% CI','97.5% CI','Cumulative sum of means','Fraction of cumulative sum','Analysed flux-microbe pair'};
mContributions = renamevars(mContributions,currNames,newNames);

% Load the analysed flux-microbe associations
filteredfmLinksSummary = readtable(filteredfmLinksSummaryPath,'VariableNamingRule','preserve');
filteredfmLinksSummary = renamevars(filteredfmLinksSummary,'Row','Microbial species');

end