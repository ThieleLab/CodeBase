% extract the results from random forests classifiers and summarize in a
% table

% define the four types of data
datatypes={
    'Strain-level abundance','Organism Abundance'
    'Reaction abundance','Reaction Abundance'
    'Reaction presence','Reaction Presence'
    'Subsystem abundance','Subsystem Abundance'
    };

defineScenarios

%% extract a table summarizing the predictive scores of each analysis

% define the four types of data
% prepare the table
scoreTable = {'Dataset',''};
cnt=2;

for i=1:size(datatypes,1)
    scoreTable{1,cnt} = datatypes{i};

    for j=1:size(datasets,1)
          % get the summary file from each analysis
        % then fill out the table with the relevant data
        scoreTable{j+1,1} = datasets{j,1};
        data=readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Subgroup_analysis' filesep 'RF_Results' filesep datasets{j,1} filesep 'summary.txt']);
        feat = find(strcmp(data(:,1),'Model Metrics Final'));
        % cut off data above the final score
        data(1:feat(end),:) = [];
        % find the score for each data type
        dt = find(strcmp(data(:,1),datatypes{i,2}));
        scoreTable{j+1,cnt} = strrep(data{dt+4,1},'Accuracy: ','');
        scoreTable{j+1,cnt} = round(str2double(scoreTable{j+1,cnt}),2,'significant');
        % get the number of reduced features
        redFeat=strsplit(data{dt+1,1},'> ');
        scoreTable{j+1,cnt} = [num2str(scoreTable{j+1,cnt}) ' (' redFeat{1,2} ')'];
        if strncmp(scoreTable{j+1,cnt},'1 ',2)
            scoreTable{j+1,cnt} = strrep(scoreTable{j+1,cnt},'1 ','>0.99 ');
        end
    end
    cnt=cnt+1;
end
writetable(cell2table(scoreTable),[rootDir filesep 'results' filesep 'microbiomes' filesep 'PredictionScoresCombined.csv'],'writeVariableNames',false)
