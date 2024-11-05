
% extract the results from random forests classifiers and summarize in a
% table

clear all
rootDir = pwd;

% define the three types of data
datatypes={
    'Internal metabolite production','internal_production'
    'Metabolite uptake and secretion','uptake_secretion'
    'Reaction presence','reaction'
    };

% define the different datasets to extract
datasets={
    'Pasolli',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'MachineLearning_Results' filesep 'Pasolli'];
    'Almeida',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'MachineLearning_Results' filesep 'Almeida'];
    'All strains',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'MachineLearning_Results' filesep 'Pasolli_Almeida_combined'];
    'Pasolli on Almeida',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'MachineLearning_Results' filesep 'Pasolli_on_Almeida'];
    'Almeida on Pasolli',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'MachineLearning_Results' filesep 'Almeida_on_Pasolli'];
    };

% define the files with the results
files={
    'feature_importance' filesep 'all_feature_importance.csv'
    'feature_importance' filesep 'final_feature_importance.csv'
    'feature_importance' filesep 'important_feature_importance.csv'
    };


%% extract a table summarizing the predictive scores of each analysis
% prepare the table
scoreTable = {
    '','Model metrics, reduced features','','','','',''
    'Dataset','Phylum','Class','Order','Family','Genus','Species'
    };
cnt=3;

for j=1:size(datasets,1)
    for i=1:size(datatypes,1)
        % first get all features in each data type analysis
        % navigate to the folders
        cd(datasets{j,2})
        % get the content of the folder
        dInfo = dir(pwd);
        folders={dInfo.name};
        folders=folders';
        folders(find(strncmp(folders,'.',1)),:)=[];
        find_folder=find(strncmp(folders,datatypes{i,2},length(datatypes{i,2})));
        if ~isempty(find_folder)
            cd(folders{find_folder})
            % read the important files from each subfolder on the taxon levels
            dInfo = dir(pwd);
            folders={dInfo.name};
            folders=folders';
            folders(find(strncmp(folders,'.',1)),:)=[];
            % get the summary file from each analysis
            % then fill out the table with the relevant data
            scoreTable{cnt,1} = [datasets{j,1} ', ' lower(datatypes{i,1})];
            for k=2:7
                getFolder=find(strcmp(folders,scoreTable{2,k}));
                data=readInputTableForPipeline([folders{getFolder} filesep 'summary.txt']);
                mmRedFeat = find(strcmp(data(:,1),'Model Metrics Final'));
                scoreTable{cnt,k} = strrep(data{mmRedFeat(1)+4,1},'Accuracy: ','');
                scoreTable{cnt,k} = round(str2double(scoreTable{cnt,k}),2,'significant');
                % get the number of reduced features
                redFeat=strsplit(data{mmRedFeat(1)+1,1},'> ');
                scoreTable{cnt,k} = [num2str(scoreTable{cnt,k}) ' (' redFeat{1,2} ')'];
                if strncmp(scoreTable{cnt,k},'1 ',2)
                    scoreTable{cnt,k} = strrep(scoreTable{cnt,k},'1 ','>0.99 ');
                end
            end
            cnt=cnt+1;
        end
        cd(rootDir)
    end
end
writetable(cell2table(scoreTable),[rootDir filesep 'results' filesep 'strains' filesep 'Table_1.csv'],'writeVariableNames',false)
