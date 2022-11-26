% generate summary tables of comparison by resource for Figure 3

% overview table of general model statistics
% load the relevant data
load([rootDir filesep 'ValidationAgainstExperimentalData' filesep 'Model_statistics.mat'])
resources=sort(fieldnames(stats));

TableStatsOverview={
    '','','','','','',''
    '# reconstructions','','','','','',''
    'Average # of reactions','','','','','',''
    'Average # of metabolites','','','','','',''
    'Average # of genes','','','','','',''
    'Average # of stoich. consistent reactions','','','','','',''
    'Average # of flux consistent reactions','','','','','',''
    'ATP yield on aerobic CM','','','','','',''
    'ATP yield on anaerobic CM','','','','','',''
    };

for i=1:length(resources)
    TableStatsOverview{1,i+1}=resources{i};
    TableStatsOverview{2,i+1}=size(stats.(resources{i}),1);
    TableStatsOverview{3,i+1}=mean(stats.(resources{i})(:,1));
    TableStatsOverview{4,i+1}=mean(stats.(resources{i})(:,2));
    TableStatsOverview{5,i+1}=mean(stats.(resources{i})(:,3));
end

load([rootDir filesep 'ValidationAgainstExperimentalData' filesep 'Stoch_Flux_Consistency.mat'])

for i=1:length(resources)
    TableStatsOverview{6,i+1}=mean(SFconsist.(resources{i})(:,2));
    TableStatsOverview{7,i+1}=mean(SFconsist.(resources{i})(:,4));
end

load([rootDir filesep 'ValidationAgainstExperimentalData' filesep 'ATP_production.mat'])

for i=1:length(resources)
    TableStatsOverview{8,i+1}=mean(atp.(resources{i})(:,1));
    TableStatsOverview{9,i+1}=mean(atp.(resources{i})(:,2));
end

cell2csv([rootDir filesep 'ValidationAgainstExperimentalData' filesep 'Results' filesep 'OverViewTableStats.csv'],TableStatsOverview);

% Table comparing results of the validation against experimental data

TableValidation={
    '','AGORA2','BiGG','CarveMe','gapseq_published','gapseq_this_study','KBase','MAGMA'
    '# models tested against NJC19','','','','','','',''
    'Accuracy against NJC19','','','','','','',''
    'Sensitivity against NJC19','','','','','','',''
    'Specificity against NJC19','','','','','','',''
    '# models tested against BacDive (metabolites)','','','','','','',''
    'Accuracy against BacDive (metabolites)','','','','','','',''
    'Sensitivity against BacDive (metabolites)','','','','','','',''
    'Specificity against BacDive (metabolites)','','','','','','',''
    '# models tested against BacDive (enzymes)','','','','','','',''
    'Accuracy against BacDive (enzymes)','','','','','','',''
    'Sensitivity against BacDive (enzymes)','','','','','','',''
    'Specificity against BacDive (enzymes)','','','','','','',''
    '# models tested against Madin','','','','','','',''
    'Sensitivity against Madin','','','','','','',''
    };

% fill out results for resources one by one
data2Load={
    [rootDir filesep 'ValidationAgainstExperimentalData' filesep 'AGORA2'],2,3
    [rootDir filesep 'ValidationAgainstExperimentalData' filesep 'BiGG'],3,4
    [rootDir filesep 'ValidationAgainstExperimentalData' filesep 'CarveMe'],4,4
    [rootDir filesep 'ValidationAgainstExperimentalData' filesep 'gapseq'],5,4
    [rootDir filesep 'ValidationAgainstExperimentalData' filesep 'gapseq_new'],6,4
    [rootDir filesep 'ValidationAgainstExperimentalData' filesep 'AGORA2'],7,2
    [rootDir filesep 'ValidationAgainstExperimentalData' filesep 'MAGMA'],8,4
    };

for i=1:size(data2Load,1)

    % define the right columns to get and save the data
tabCol=data2Load{i,2};
resCol=data2Load{i,3};

% load the datasets one by one
load([data2Load{i,1} filesep 'Results_NJC19.mat'])
TableValidation{2,tabCol}=Results{9,resCol};
TableValidation{3,tabCol}=Results{8,resCol};
TableValidation{4,tabCol}=Results{6,resCol};
TableValidation{5,tabCol}=Results{7,resCol};

load([data2Load{i,1} filesep 'Results_BacDive_Metabolites.mat'])
TableValidation{6,tabCol}=Results{9,resCol};
TableValidation{7,tabCol}=Results{8,resCol};
TableValidation{8,tabCol}=Results{6,resCol};
TableValidation{9,tabCol}=Results{7,resCol};

load([data2Load{i,1} filesep 'Results_BacDive_Enzymes.mat'])
TableValidation{10,tabCol}=Results{9,resCol};
TableValidation{11,tabCol}=Results{8,resCol};
TableValidation{12,tabCol}=Results{6,resCol};
TableValidation{13,tabCol}=Results{7,resCol};

load([data2Load{i,1} filesep 'Results_Madin.mat'])
TableValidation{14,tabCol}=Results{5,resCol};
TableValidation{15,tabCol}=Results{4,resCol};
end

cell2csv([rootDir filesep 'ValidationAgainstExperimentalData' filesep 'Results' filesep 'OverViewTableValidation.csv'],TableValidation);

