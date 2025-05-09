function stats = makeMetadataTable(metadata)
% Personal function for the Parkinson project by Tim Hensen.
% This function creates a summary statistics table for the manuscript. 
% 

% Load the metadataFiltered and preallocate an empty table and then populate the table with number of samples

% Load metadataFiltered in memory
% metadata = readtable(metadataPath);

% Remove variables that will not be used in the metadataFiltered table
varsToRemove = {'ID','collection_method','total_sequences','Hispanic_or_Latino','Race','Jewish_ancestry','Bristol_stool_chart'};
metadataFiltered = metadata;
metadataFiltered(:,varsToRemove) = [];

% Process dietary information and make it binary. Only look at daily
% intake. 

% Get variables with dietary informations
metadataVars = metadataFiltered.Properties.VariableNames;
dietVars = metadataVars(contains(metadataVars,'How_often_do_you_eat_') & ~contains(metadataVars,'YOGURT')); % Process yoghurt separately
dietInfo = metadataFiltered(:,dietVars);

% Find individuals that consume a given food category daily
dietInfoArray = table2array(dietInfo);
dietInfoArray(matches(dietInfoArray,{'Few times a week','Few times a month','Less than once a month or never'})) = {'N'};
dietInfoArray(matches(dietInfoArray,'At least once a day')) = {'Y'};
metadataFiltered{:,dietVars} = dietInfoArray;

% Update variable names
dailyDietVars = {'Grains, daily','Poultry, beef, pork, seafood, or eggs, daily','Fruits, daily','Nuts, daily'};
metadataFiltered.Properties.VariableNames(matches(metadataVars,dietVars)) = dailyDietVars;

% Next, process Yoghurt intake 
yogurtIntake = metadataFiltered.How_often_do_you_eat_YOGURT;
yogurtIntake(matches(yogurtIntake,{'At least once a day','Few times a month','Less than once a month or never'})) = {'N'};
yogurtIntake(matches(yogurtIntake,'Few times a week')) = {'Y'};
metadataFiltered.How_often_do_you_eat_YOGURT = yogurtIntake;
metadataFiltered = renamevars(metadataFiltered,"How_often_do_you_eat_YOGURT","Yogurt, a few times a week");

% Next, check which variables have only one category and remove these
% variables from the table
% Find numerical variables and not test them
notAnumVars = ~table2array(varfun(@isnumeric,metadataFiltered));

% Also identify variables wit only N or missing information and remove them
metadataFilteredNoNumVars = metadataFiltered(:,notAnumVars);
varsToRemove = table2array(varfun(@(x) ~any(matches(x,{'Y','PD','male'})), metadataFilteredNoNumVars));
size(metadataFiltered)
varNames = metadataFilteredNoNumVars.Properties.VariableNames;
varNamesOfRemovedVars = varNames(varsToRemove);
metadataFiltered(:,varNamesOfRemovedVars)=[];
size(metadataFiltered)
%

% Preallocate table with variables of interest
variables = metadataFiltered.Properties.VariableNames;
stats = array2table(nan(length(variables),6));
stats.Properties.VariableNames = {'PD_n','PD_stat','Control_n','Control_stat','p','OR'};
stats.Properties.RowNames = variables';%[{'Samples'} variables]';

% Populate the table with the number of samples using the local sampleCount
% function
for i=1:length(variables)
    stats{stats.Properties.RowNames{i},{'PD_n','Control_n'}} = sampleCount(metadataFiltered,variables{i});
end
%
% Obtain the following summary statistics for the variables of interest
% mean + SD for age and BMI

% Make the indices strings
stats = convertvars(stats,{'PD_stat','Control_stat','p','OR'},'string');

[groups,groupNames] = findgroups(metadataFiltered.Case_status);

% Find the mean age for pd and controls
meanSd = @(x) string(round([mean(x,'omitnan') std(x,'omitnan')],2));
ageMeanSd = splitapply(meanSd, metadataFiltered.Age_at_collection,groups) ;
bmiMeanSd = splitapply(meanSd, metadataFiltered.BMI,groups);

% Convert results for latex
convertFun = @(x,y) append(x(matches(groupNames,y),1)," $\pm$ ", x(matches(groupNames,y),2));

% Add results to table
stats{'Age_at_collection','PD_stat'} = convertFun(ageMeanSd,'PD');
stats{'Age_at_collection','Control_stat'} = convertFun(ageMeanSd,'Control');
stats{'BMI','PD_stat'} = convertFun(bmiMeanSd,'PD');
stats{'BMI','Control_stat'} = convertFun(bmiMeanSd,'Control');

% Find PD and Control
pd = matches(metadataFiltered.Case_status,'PD');
cn = matches(metadataFiltered.Case_status,'Control');

% Add p-values for numerical data using the ranksum wilcoxon test
stats{'Age_at_collection','p'} = string(sprintf('%.2e', ranksum(metadataFiltered.('Age_at_collection')(pd),metadataFiltered.('Age_at_collection')(cn))));
stats{'BMI','p'} = string(sprintf('%.2e', ranksum(metadataFiltered.('BMI')(pd),metadataFiltered.('BMI')(cn))));

%
% Find the number of samples and percentage of true values
% Add count and percentage for binary data

% For sex, convert male/female to Y/N
metadataFiltered.Sex(matches(metadataFiltered.Sex,'male'))={'Y'};
metadataFiltered.Sex(matches(metadataFiltered.Sex,'female'))={'N'};

% Get the binary variables, including sex
binaryVars = variables(~matches(variables,{'Age_at_collection','BMI'}))';

% Find the number of samples with Y and N
sampleStats = @(data,value) sum(matches(data,value))';
roundConvert = @(x) string(round(x,2));

% Find the number of samples with the condition for PD patients
pdData = metadataFiltered{pd,binaryVars};
pdNums = sampleStats(pdData,'Y');
pdNums(:,2) = pdNums ./ (pdNums + sampleStats(pdData,'N')) * 100;
pdNums = roundConvert(pdNums);

% Find the number of samples with the condition for Controls
cnData = metadataFiltered{cn,binaryVars};
cnNums = sampleStats(cnData,'Y');
cnNums(:,2) = cnNums ./ (cnNums + sampleStats(cnData,'N')) * 100;
cnNums = roundConvert(cnNums);

% Convert the data for latex table and append to table
stats{binaryVars,'PD_stat'} = append(pdNums(:,1), " ( ",pdNums(:,2)," )");
stats{binaryVars,'Control_stat'} = append(cnNums(:,1), " ( ",cnNums(:,2)," )");
%
% Obtain p-values for categorical variables
%format shortEng
for i=1:length(binaryVars)

    % Remove samples without information
    metadataTemp = metadataFiltered;

    % Remove samples with no information on the binary variable of interest
    metadataTemp(cellfun(@isempty, metadataFiltered.(binaryVars{i})),:) = [];
    metadataTemp = metadataTemp(:,{'Case_status',binaryVars{i}});
    
    % Crosstabulate
    cross = crosstab(metadataTemp.Case_status,metadataTemp.(binaryVars{i}));

    % Calculate fisher p-value and add to table
    [~,p,ORstats] = fishertest(cross);

    % Format p-value
    stats{binaryVars{i},'p'} = string(sprintf('%.2e', p));

    % Add odds ratio and 95% ci to table
    ORstats = structfun(@(x) string(round(x,2)),ORstats,'UniformOutput',false);

    % Add OR and ci to table
    stats{binaryVars{i},'OR'} = strcat(ORstats.OddsRatio," [",ORstats.ConfidenceInterval(1),"—",ORstats.ConfidenceInterval(2),"]");
end
%

% Make all columns a string
stats = convertvars(stats,stats.Properties.VariableNames,"string");

% Set all OR results that include infinities to ""
stats.OR(contains(stats.OR,"Inf"))=missing;

% Change table for latex compiling

% Variables
stats.Properties.VariableNames = regexprep(stats.Properties.VariableNames,'_',' ');
stats.Properties.VariableNames = regexprep(stats.Properties.VariableNames,'p','\\textit{P}');
stats.Properties.VariableNames = regexprep(stats.Properties.VariableNames,'stat','summary statistic');
stats.Properties.VariableNames = regexprep(stats.Properties.VariableNames,' n',' \\textit{N}');
stats.Properties.VariableNames = regexprep(stats.Properties.VariableNames,'OR','OR [95\\% CI]');

% Process rows

% Remove case_status row
stats('Case_status',:) = [];

% Order rows so that sex, age, and bmi are at the top
newRowOrder = stats.Properties.RowNames(~matches(stats.Properties.RowNames,{'Sex','Age_at_collection','BMI'}));
newRowOrder = [{'Sex';'Age_at_collection';'BMI'} ; newRowOrder];
[~,~,ib] = intersect(newRowOrder,stats.Properties.RowNames,'stable');
stats = stats(ib,:);

% Now process the row names

% Improve row legibility
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Age_at_collection','Mean age \\pm SD');
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'BMI','Mean BMI \\pm SD');
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Sex','Male samples (\\% male)');

% Remove day of stool collection substring
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Day_of_stool_collection_','');

% Process weight gain/loss names
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Loss_','Lost_');

% Process Lifestyle questions
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Do_you_smoke','Smoking');
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Do_you_drink_caffeinated_beverages','Caffeine consumption');
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Do_you_drink_alcohol','Alcohol consumption');

% Process drug and medication names 
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'_drugs',' drugs');
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'_med',' medications');

 % Remove underscores
stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'_',' ');


% stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Depression anxiety mood med','Mood medication users (\\% users)');
% stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Diabetes med','Diabetes medication users (\\% users)');
% stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Sleep aid','Sleep aid medication users (\\% users)');
% stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Cholesterol med','Cholesterol medication users (\\% users)');
% stats.Properties.RowNames = regexprep(stats.Properties.RowNames,'Blood pressure med','Blood pressure medication users (\\% users)');

% Save table
% metadataStatsPath = regexprep(metadataPath, '[\\/][^\\/]*$', '');
% metadataStatsPath = [char(metadataStatsPath) filesep 'metadataTable.xlsx'];
% 
% writetable(stats,metadataStatsPath,'WriteRowNames',true)
end

function N = sampleCount(metadataFiltered,variable)
% Get the number of samples
variable = string(variable);

% Get PD and controls
pd = matches(metadataFiltered.Case_status,'PD');
cn = matches(metadataFiltered.Case_status,'Control');

if isa(metadataFiltered.(variable),'double')
    N = [sum(~isnan(metadataFiltered.(variable)(pd))) sum(~isnan(metadataFiltered.(variable)(cn)))];
elseif matches(variable,"Case_status")
    N = [sum(pd) sum(cn)];
elseif contains(variable,"Sex")
    N = [sum(matches(metadataFiltered.(variable)(pd),{'male','female'})) sum(matches(metadataFiltered.(variable)(cn),{'male','female'}))];
else
    N = [sum(matches(metadataFiltered.(variable)(pd),{'Y','N'})) sum(matches(metadataFiltered.(variable)(cn),{'Y','N'}))];
end

end
