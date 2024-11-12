function [genOverview] = statAnalysisSubgroups(subgroup, stratVar)

%%statistical analysis on 240K subgroups

%% Read metadata file
% Import the metadata of the subgroup
metadata=readtable(append(subgroup,'_samples.csv'),'ReadVariableNames',true);
try
    % Remove any row that is already called sample, not important for
    % analysis
    metadata = removevars(metadata, {'Sample'});
catch
end

% Rename first column to Sample to ensure it can be combined with other
% tables later on
metadata.Properties.VariableNames{1} = 'Sample';

% Ensure that column 1 is the correct data type
metadata = convertvars(metadata, metadata.Properties.VariableNames(1),'string');
% Obtain amount of samples for each unique group in the stratification
% variable
genOverview = groupsummary(metadata, stratVar);

% Obtain the sample names
sampleNames = metadata.Sample;

%% Read model properties

% Read model properties files
modelProperties = readtable(append('Model_sizes_',subgroup,'.csv'),'Readrownames',1);

% Transpose table
modelProperties = rows2vars(modelProperties);

% Set first column name to Sample
modelProperties.Properties.VariableNames{1} = 'Sample';

% Ensure that column 1 is the correct data type
modelProperties = convertvars(modelProperties, modelProperties.Properties.VariableNames(1),'string');

%% Read relative abundance microbial species
% Read the microbial relative abundances file
microbeAbnTmp= readInputTableForPipeline(append('Organism_abundance_',subgroup,'.csv'));

% Transpose and convert to table
microbeAbn = cell2table(microbeAbnTmp(:,2:end)', 'VariableNames',microbeAbnTmp(:,1));

% Set first column name to Sample
microbeAbn.Properties.VariableNames{1} = 'Sample';

% Ensure correct data types in table
microbeAbn = convertvars(microbeAbn, microbeAbn.Properties.VariableNames(1),'string');
microbeAbn = convertvars(microbeAbn, microbeAbn.Properties.VariableNames(2:end),'double');

%% Read the reaction relative abundance file

% Read the reaction abundance table
rxnAbnTmp = readInputTableForPipeline(append('Reaction_abundance_',subgroup,'.csv'));

% Transpose and convert to table
rxnAbn = cell2table(rxnAbnTmp(:,2:end)', 'VariableNames',rxnAbnTmp(:,1));

% Set first column name to Sample
rxnAbn.Properties.VariableNames{1} = 'Sample';

% Ensure correct data types in table
rxnAbn= convertvars(rxnAbn, rxnAbn.Properties.VariableNames(1),'string');
rxnAbn = convertvars(rxnAbn, rxnAbn.Properties.VariableNames(2:end),'double');

%% Read the reaction presence file
% Read the reaction presence table
rxnPrsnTmp = readInputTableForPipeline(append('Reaction_presence_',subgroup,'.csv'));

% Transpose and convert to table
rxnPrsn = cell2table(rxnPrsnTmp(:,2:end)', 'VariableNames',rxnPrsnTmp(:,1));

% Set first column name to Sample
rxnPrsn.Properties.VariableNames{1} = 'Sample';

% Ensure correct data types in table
rxnPrsn = convertvars(rxnPrsn, rxnPrsn.Properties.VariableNames(1),'string');
rxnPrsn = convertvars(rxnPrsn, rxnPrsn.Properties.VariableNames(2:end),'double');

%% Read the subsystem relative abundance file

% Import subsystem relative abundances for each sample
subsysAbn = readtable(append('Subsystem_abundance_',subgroup,'.csv'),'Readrownames',1);

% Transpose table
subsysAbn=rows2vars(subsysAbn);

% Set first column name to Sample
subsysAbn.Properties.VariableNames{1} = 'Sample';

% Ensure correct data types in table
subsysAbn = convertvars(subsysAbn, subsysAbn.Properties.VariableNames(1),'string');

%% Read in the simulation results if performed

% Check if the Objective.txt file exists
if isfile('Objectives_AED.txt')
    % Import reaction relative abundance for each sample
    % readInputTableForPipeline used instead of readtable to ensure correct
    % reaction names
    fluxTmp = readInputTableForPipeline(append('Objectives_AED.txt'));
    
    % Transpose and convert to table
    fluxes = cell2table(fluxTmp(:,2:end)', 'VariableNames',fluxTmp(:,1));
    
    % Set first column name to Sample
    fluxes.Properties.VariableNames{1} = 'Sample';
    
    % Ensure correct data types in table
    fluxes = convertvars(fluxes, fluxes.Properties.VariableNames(1),'string');
    fluxes = convertvars(fluxes, fluxes.Properties.VariableNames(2:end),'double');
    
    % Remove microbiota_model_samp_ from the sample names in the flux table
    % to ensure they can be matched with the metadata
    fluxes.Sample = strrep(fluxes.Sample, 'microbiota_model_samp_', '');
    
end
%% Load the VMH database to obtain the various subsystems for reactions
rxnSubDatabase = loadVMHDatabase;
% Extract the reaction table
rxnSubDatabase = rxnSubDatabase.reactions;

% Extract the reaction ID and the two subsystem columns
rxnSubDatabase = rxnSubDatabase(2:end,[1,11,12]);


%% Clear not required variables
clear rxnPrsnTmp rxnAbnTmp microbeAbnTmp

%% Check for number of stratification variables
% Extract the column with the stratification variable form the metadata
stratArray = metadata.(stratVar);

% Determine if there are 2 or more unique groups for stratification
if size(unique(stratArray),1) > 2
    anovaBool = true;
else
    anovaBool = false;
end

%% generate index vectors for grouping variables

if anovaBool
    % For anova, generate a array with each sample ID corresponding to the
    % stratification name
    anovaArray = cat(2,sampleNames, stratArray);
    if ~istable(anovaArray)
        anovaArray = array2table(anovaArray);
        anovaArray.Properties.VariableNames = {'Sample', 'stratVar'};
    else
        anovaArray.Properties.VariableNames(2) = {'stratVar'};
    end
    
    % To avoid errors
    stratCatNames = {};
    catStratVarHeader = {};
else
    % For only two groups, change the stratifcation to 0/1 instead of
    % strings
    catStratVarHeader = ['cat', stratVar];
    
    % Convert variables into categories and save in table
    [stratCat, stratCatNames] = grp2idx(metadata.(stratVar));
    metadata.(catStratVarHeader)= stratCat;
    
    % Convert the header name into cell
    catStratVarHeader = {catStratVarHeader};
    
    % Covert 1/2 variables to 0/1
    metadata.(string(catStratVarHeader{1}))=metadata.(string(catStratVarHeader{1}))-1;
end

%% remove reactions, microbe and subsysem abundance variables with SD=0

% find the width of the table
n=size(rxnAbn,2);

% stop at n = 1 as that column contains strings. The investigation starts
% at the end of the table and moves backwards
while n>=2
    
    % calculate the standard deviation and if it is 0 remove the column
    % from the table
    if std(rxnAbn.(n))==0
        rxnAbn.(n)=[];
    end
    
    % move to the preceding column
    n=n-1;
end

% find the width of the table
n=size(subsysAbn,2);

% stop at n = 1 as that column contains strings. The investigation starts
% at the end of the table and moves backwards
while n>=2
    
    % calculate the standard deviation and if it is 0 remove the column
    % from the table
    if std(subsysAbn.(n))==0
        subsysAbn.(n)=[];
    end
    
    % move to the preceding column
    n=n-1;
end

% find the width of the table
n=size(microbeAbn,2);

% stop at n = 1 as that column contains strings. The investigation starts
% at the end of the table and moves backwards
while n>=2
    
    % calculate the standard deviation and if it is 0 remove the column
    % from the table
    if std(microbeAbn.(n))==0
        microbeAbn.(n)=[];
    end
    
    % move to the preceding column
    n=n-1;
end

%% remove categorical reactions present in less than 5% or more than 95% of the models

% find the width of the table
n=size(rxnPrsn,2);

% stop at n = 1 as that column contains strings. The investigation starts
% at the end of the table and moves backwards
while n>=2
    % Calculate the mean for each column and check if the mean is lower
    % than 0.05 (5%) or higher than 0.95 (95%)
    if mean(rxnPrsn.(n))<0.05||mean(rxnPrsn.(n))>0.95
        rxnPrsn.(n)=[];
    end
    % move to the preceding column
    n=n-1;
end
%% Calculate diversity measures

% Initialise a new table with the first column being the sample names
diversity=microbeAbn(:,1);

% Initialise a new column for alpha-diversity
diversity.alpha(1,1)=0;

% Calculate the alpha-diversity
for i=2:size(microbeAbn,2)
    diversity.alpha=diversity.alpha+(microbeAbn.(i)>0);
end

% Initialise a new column for beta-diversity
diversity.betastand(1,1)=0;

% Calculate the beta-diversity standardized for the number of samples
diversity.betastand = (size(microbeAbn,2)./diversity.alpha)./size(microbeAbn,1);

% Initialise columns for Shannon (shannon), Shannon Evennes (SEI), Simpson diversity
% index (SDI)
diversity.shannon(1,1)=0;
diversity.SEI(1,1)=0;

% This column is an inverse of the SDI and is used to calculate the SDI
diversity.simpson(1,1)=0;
diversity.SDI(1,1)=0;

% Calculate simpson and Shannon Diversity
for i=2:size(microbeAbn,2)
    diversity.simpson = diversity.simpson+(microbeAbn.(i).*microbeAbn.(i));
    
    singleSpeciesShannonCoeff = microbeAbn.(i).*log(microbeAbn.(i));
    singleSpeciesShannonCoeff(isnan(singleSpeciesShannonCoeff)) = 0;
    diversity.shannon=diversity.shannon-singleSpeciesShannonCoeff;
end

% Calculate the SEI and SDI
diversity.SEI=diversity.shannon./log(diversity.alpha);
diversity.SDI=1-diversity.simpson;

%% statistical analysis on models properties and diversity measures

% Check if the Readcount variable is present in the table. If yes include
% in analysis if no, leave out of analysis
if sum(strcmp('Readcount', metadata.Properties.VariableNames)) == 1
    % Create headers for the table to store data in
    tableHeaders = [{'ReadCount'} diversity.Properties.VariableNames(:,2:end) modelProperties.Properties.VariableNames(:,2:end)];
else
    % Create headers for the table to store data in
    tableHeaders = [diversity.Properties.VariableNames(:,2:end) modelProperties.Properties.VariableNames(:,2:end)];
end

% Check if there are only 2 or more stratification groups
if anovaBool
    modelStatTable = join(anovaArray,modelProperties);
    modelStatTable = join(modelStatTable, diversity);
    modelStatTable = join(modelStatTable, metadata(:,[1,2]));
    
    % Create names for table headers
    firstTabVarNames = {'mod','Skewness','Kurtosis'};
else
    modelStatTable=join(metadata,modelProperties);
    modelStatTable=join(modelStatTable,diversity);
end

% Generate filename
filename=[strcat(subgroup,'_stat_mod.csv')];
% run statistics

continuousValuesStatistics(modelStatTable, tableHeaders, 'mod', stratCatNames, filename, rxnSubDatabase, anovaBool, catStratVarHeader)

%% statistical analysis on relative abundances of microbes
% Obtain all microbess (the table column names)

microbeNames=microbeAbn.Properties.VariableNames(:,2:end);
% Create the table to investigate depending if there are 2 or more
% stratification groups
if anovaBool
    modelStatTable = join(anovaArray,microbeAbn);
else
    modelStatTable=join(metadata,microbeAbn);
end

% Generate filename
filename=[strcat(subgroup,'_stat_micr.csv')];
% run statistics

continuousValuesStatistics(modelStatTable, microbeNames, 'microbes', stratCatNames, filename, rxnSubDatabase, anovaBool, catStratVarHeader)
%% statistical analysis on relative abundances of reactions
% Obtain all reaction names (the table column names)
rxn=rxnAbn.Properties.VariableNames(:,2:end);
% Create the table to investigate depending if there are 2 or more
% stratification groups
if anovaBool
    modelStatTable = join(anovaArray,rxnAbn);
else
    modelStatTable=join(metadata,rxnAbn);
end
% Generate filename

filename=[strcat(subgroup,'_stat_rxn.csv')];
% run statistics

continuousValuesStatistics(modelStatTable, rxn, 'reactions', stratCatNames, filename, rxnSubDatabase, anovaBool, catStratVarHeader)
%% statistical analysis on relative abundances of subsystems
% Obtain all subsystems (column headers)
subystem=subsysAbn.Properties.VariableNames(:,2:end);
% Create the table to investigate depending if there are 2 or more
% stratification groups
if anovaBool
    modelStatTable = join(anovaArray,subsysAbn);
else
    modelStatTable=join(metadata,subsysAbn);
end
% Generate filename

filename=[strcat(subgroup,'_stat_subsys.csv')];
% run statistics

continuousValuesStatistics(modelStatTable, subystem, 'subsystems', stratCatNames, filename, rxnSubDatabase, anovaBool, catStratVarHeader)
%% statistical analysis on microbial fluxes
% Check if the file with flux values exists
if isfile('Objectives_AED.txt')
    % Obtain all flux names (column names)
    fluxNames=fluxes.Properties.VariableNames(:,2:end);
    
    % Create the table to investigate depending if there are 2 or more
    % stratification groups
    if anovaBool
        modelStatTable = join(anovaArray,fluxes);
    else
        modelStatTable=join(metadata,fluxes);
    end
    % generate filename
    filename=[strcat(subgroup,'_stat_flux.csv')];
    
    % run statistics
    continuousValuesStatistics(modelStatTable, fluxNames, 'excretionfluxes', stratCatNames, filename, rxnSubDatabase, anovaBool, catStratVarHeader)
end

%% statistical analysis on absence presence of reactions

% Obtain array with 2:end column headers
rxnID=rxnPrsn.Properties.VariableNames(:,2:end);

if anovaBool
    modelStatTable = join(anovaArray,rxnPrsn);
    % Calculate for each reaction the total amount of absent and present
    % for each subgroup and the p-value for Pearson Chi-squared
    for i = 1:length(rxnID)
        [tbl, chi2,p,labels] = crosstab(modelStatTable.stratVar, modelStatTable.(rxnID{i}));
        presenceArray = {'x'};
        for j = 1:size(labels,1)
            presenceArray(1,end+1:end+2) = {tbl(j,1),tbl(j,2)};
        end
        
        % Find the reaction ID in the database
        f = find(contains(rxnSubDatabase(:,1), rxnID{i}));
        
        % If more option are available use the strcmp function for an
        % exact match
        if length(f) > 1
            f = find(strcmp(rxnSubDatabase(:,1), rxnID{i}));
            subsystems(1,1:2) = rxnSubDatabase(f, 2:3);
            
            % If no matches are found, remove the last character from the
            % string and search for an exact match (Due to a MATLAB
            % formatting error that adds on characters to some strings)
        elseif isempty(f)
            rxn = rxnID{i};
            rxn = rxn(1:end-1);
            f = find(strcmp(rxnSubDatabase(:,1), rxn));
            subsystems(1,1:2) = rxnSubDatabase(f, 2:3);
            % Change rxn ID to correct form
            rxnID{i} = rxn;
        else
            subsystems(1,1:2) = rxnSubDatabase(f, 2:3);
        end
        % Store data in an array
        storeArray(i,:) = {subsystems{1,:},presenceArray{1,2:end}, p};
    end
    % Obtain table column names
    firstTblVarNames = {'x'};
    for k = 1:size(labels,1)
        firstTblVarNames(1,end+1:end+2) = {strcat('Absent', labels{k,1}), strcat('Present', labels{k,1})};
    end
    tableVarNames = ['reactions', 'subsystems', 'subsystems_general',firstTblVarNames(2:end), 'p_PearshonChi2'];
    % Create table
    rxnPrsnStatTable = table(rxnID');
    rxnPrsnStatTable(:,2:(size(storeArray,2)+1)) = storeArray;
    rxnPrsnStatTable.Properties.VariableNames = tableVarNames;
    rxnPrsnStatTable=sortrows(rxnPrsnStatTable,size(rxnPrsnStatTable,2));
    
    totalTests=size(rxnPrsnStatTable,1);
    rxnPrsnStatTable.FDRCorrection=totalTests*table2array(rxnPrsnStatTable(:,end));
else
    modelStatTable=join(metadata,rxnPrsn);
    for i=1:length(catStratVarHeader)
        for j=1:length(rxnID)
            % Calculate for each reaction the total amount of absent and present
            % for each subgroup and the p-value for Pearson Chi-squared
            [tbl,chi2,p,labels] = crosstab(modelStatTable.(string(catStratVarHeader{i})),modelStatTable.(string(rxnID{j})));
            var0ifrxnabsent{j}=tbl(1,1);
            var1ifrxnabsent{j}=tbl(2,1);
            var0ifrxnpresent{j}=tbl(1,2);
            var1ifrxnpresent{j}=tbl(2,2);
            % Choose here to do chi2 or fisher test for 240k study
            %             p_chi{j}=p;
            
            [~, pval, ~] = fishertest(tbl);
            p_fisher{j}=pval;
            
            % Find the reaction ID in the database
            f = find(contains(rxnSubDatabase(:,1), rxnID{j}));
            
            % If more option are available use the strcmp function for an
            % exact match
            if length(f) > 1
                f = find(strcmp(rxnSubDatabase(:,1), rxnID{j}));
                subsystems(j,1:2) = rxnSubDatabase(f, 2:3);
                
                % If no matches are found, remove the last character from the
                % string and search for an exact match (Due to a MATLAB
                % formatting error that adds on characters to some strings)
            elseif isempty(f);
                rxn = rxnID{j};
                rxn = rxn(1:end-1);
                f = find(strcmp(rxnSubDatabase(:,1), rxn));
                subsystems(j,1:2) = rxnSubDatabase(f, 2:3);
                % Change rxn ID to correct form
                rxnID{j} = rxn;
            else
                subsystems(j,1:2) = rxnSubDatabase(f, 2:3);
            end
        end
        
        % Create table
        rxnPrsnStatTable=table(rxnID', subsystems(:,1), subsystems(:,2),var0ifrxnabsent', var1ifrxnabsent', var0ifrxnpresent', var1ifrxnpresent', p_fisher', 'VariableNames', {'reactions', 'Subsystem','Subystem_general',strcat('absent',stratCatNames{1}), strcat('absent',stratCatNames{2}), strcat('present', stratCatNames{1}), strcat('present', stratCatNames{2}), 'p_Fisherex'});
        % TODO: Fix this so its the last row which should be p-value
        rxnPrsnStatTable=sortrows(rxnPrsnStatTable,size(rxnPrsnStatTable,2));
        
        % Adjust p-value according to Bonferonni method
        totalTests=size(rxnPrsnStatTable,1);
        rxnPrsnStatTable.FDRCorrection=totalTests*cell2mat(table2array(rxnPrsnStatTable(:,end)));
    end
end

% Save table to .csv file
filename=[strcat(subgroup,'_stat_rxnpr.csv')];
writetable(rxnPrsnStatTable,filename);
end

function continuousValuesStatistics(statisticsTable, column2investigate, dataOrigin, stratCatNames, filename, rxnSubDatabase, anovaBool, catStratVarHeader)

if anovaBool
    for i = 1:length(column2investigate)
        % Obtain the mean and standard deviations for each column to
        % investigate
        meanStd = groupsummary(statisticsTable, 'stratVar', {'mean', 'std'}, column2investigate(i));
        
        % Initialise storage variables
        meansstd = {'x'};
        partialTabVarNames = {'x'};
        for j = 1:size(meanStd,1)
            % Store the means and standard deviations in a readable way
            meansstd(1,end+1:end+2) = table2cell(meanStd(j,[3,4]));
            partialTabVarNames(1,end+1:end+2) = {strcat('mean ', meanStd{j,1}),strcat('std ', meanStd{j,1})};
        end
        
        % Calculate skewness, kurtosis and p-value column to investigate
        sk = skewness(statisticsTable.(cell2mat(column2investigate(i))));
        kur = kurtosis(statisticsTable.(cell2mat(column2investigate(i))));
        p = anova1(statisticsTable.(cell2mat(column2investigate(i))), statisticsTable.stratVar,'off');
        
        % Find the appropiate subsystems if reactions are investigated
        if strcmp(dataOrigin, 'reactions')
            % Find the reaction ID in the database
            f = find(contains(rxnSubDatabase(:,1), column2investigate{i}));
            if length(f) > 1
                
                f = find(strcmp(rxnSubDatabase(:,1), column2investigate{i}));
                subsystems(1,1:2) = rxnSubDatabase(f, 2:3);
                
                % If no matches are found, remove the last character from the
                % string and search for an exact match (Due to a MATLAB
                % formatting error that adds on characters to some strings)
            elseif isempty(f)
                rxnT = column2investigate{i};
                rxnT = rxnT(1:end-1);
                f = find(strcmp(rxnSubDatabase(:,1), rxnT));
                subsystems(1,1:2) = rxnSubDatabase(f, 2:3);
                % Change rxn ID to correct form
                column2investigate{i} = rxnT;
            else
                subsystems(1,1:2) = rxnSubDatabase(f, 2:3);
            end
            
            % Store results
            storeArray(i,:) = {subsystems{1,:}, sk, kur, meansstd{1,2:end},p};
            % Create final table headers
            tableVarNames = [{dataOrigin,'subsytem','subsystem_general','Skewness','Kurtosis'}, partialTabVarNames(1,2:end), {'pANOVA'};];
        else
            
            % Store the results
            storeArray(i,:) = {sk, kur, meansstd{1,2:end},p};
            % Create final table headers
            tableVarNames = [{dataOrigin ,'Skewness','Kurtosis'}, partialTabVarNames(1,2:end), {'pANOVA'};];
        end
    end
    
    %Create final table
    finalTable = table(column2investigate');
    finalTable(:,2:(size(storeArray,2)+1)) = storeArray;
    finalTable.Properties.VariableNames = [tableVarNames{:}];
    finalTable=sortrows(finalTable,size(finalTable,2));
else
    for i=1:length(column2investigate)
        % Obtain means, standard deviations, skewness and kurtosis
        meangroup0{i}=mean(statisticsTable.(string(column2investigate{i}))(statisticsTable.(string(catStratVarHeader))==0), 'omitnan');
        sdgroup0{i}=std(statisticsTable.(string(column2investigate{i}))(statisticsTable.(string(catStratVarHeader))==0), 'omitnan');
        meangroup1{i}=mean(statisticsTable.(string(column2investigate{i}))(statisticsTable.(string(catStratVarHeader))==1), 'omitnan');
        sdgroup1{i}=std(statisticsTable.(string(column2investigate{i}))(statisticsTable.(string(catStratVarHeader))==1), 'omitnan');
        sk{i}=skewness(statisticsTable.(string(column2investigate{i})));
        kur{i}=kurtosis(statisticsTable.(string(column2investigate{i})));
        
        % Calculate the p-value
        [~,p]=ttest2(statisticsTable.(string(column2investigate{i}))(statisticsTable.(string(catStratVarHeader))==0), statisticsTable.(string(column2investigate{i}))(statisticsTable.(string(catStratVarHeader))==1));
        p_ttest{i}=p;
        
        % Find the appropiate subsystems if reactions are investigated
        if strcmp(dataOrigin, 'reactions')
            % Find the reaction ID in the database
            f = find(contains(rxnSubDatabase(:,1), column2investigate{i}));
            
            if length(f) > 1
                f = find(strcmp(rxnSubDatabase(:,1), column2investigate{i}));
                subsystems(i,1:2) = rxnSubDatabase(f, 2:3);
                
                % If no matches are found, remove the last character from the
                % string and search for an exact match (Due to a MATLAB
                % formatting error that adds on characters to some strings)
            elseif isempty(f)
                rxnT = column2investigate{i};
                rxnT = rxnT(1:end-1);
                f = find(strcmp(rxnSubDatabase(:,1), rxnT));
                subsystems(i,1:2) = rxnSubDatabase(f, 2:3);
                % Change rxn ID to correct form
                column2investigate{i} = rxnT;
            else
                subsystems(i,1:2) = rxnSubDatabase(f, 2:3);
            end
        end       
    end
    % Create table based if the subsytems for reactions need to be added or
    % not
    if strcmp(dataOrigin, 'reactions')
        
        finalTable=table(column2investigate', subsystems(:,1), subsystems(:,2), sk', kur', meangroup0', sdgroup0', meangroup1', sdgroup1', p_ttest', 'VariableNames', {'reactions', 'Subsystem','Subystem_general','skewness', 'kurtosis',strcat('mean',stratCatNames{1}),  strcat('std',stratCatNames{1}),  strcat('mean',stratCatNames{2}),  strcat('std',stratCatNames{2}), 'p_ttest'});
    else
        finalTable=table(column2investigate', sk', kur', meangroup0', sdgroup0', meangroup1', sdgroup1', p_ttest', 'VariableNames', {dataOrigin,'skewness', 'kurtosis',strcat('mean',stratCatNames{1}),  strcat('std',stratCatNames{1}),  strcat('mean',stratCatNames{2}),  strcat('std',stratCatNames{2}), 'p_ttest'});
    end
    finalTable=sortrows(finalTable,size(finalTable,2));
end

% Adjust p-value according to Bonferonni method
totalTestsDone=size(finalTable,1);
try
   finalTable.FDRCorrection=totalTestsDone*cell2mat(table2array(finalTable(:,end)));
catch
    finalTable.FDRCorrection=totalTestsDone*table2array(finalTable(:,end));
end

% Save table as .csv file
writetable(finalTable,filename);

end
