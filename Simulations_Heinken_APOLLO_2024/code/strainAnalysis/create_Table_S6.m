
% extract classifying features for Table S6

clear all
rootDir = pwd;

table = {};
cnt=1;

% define the three types of data
datatypes={
    'Internal metabolite production','internal_production'
    'Metabolite uptake and secretion','uptake_secretion'
    'Reaction presence','reaction'
    };

% define the different datasets to extract
datasets={
    'Pasolli',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'RandomForests_Results' filesep 'Pasolli'];
    'Almeida',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'RandomForests_Results' filesep 'Almeida'];
    'All strains',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'RandomForests_Results' filesep 'Pasolli_Almeida_combined'];
    'Pasolli on Almeida',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'RandomForests_Results' filesep 'Pasolli_on_Almeida'];
    'Almeida on Pasolli',[rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'RandomForests_Results' filesep 'Almeida_on_Pasolli'];
    };

for i=1:size(datatypes,1)
    % first get all features in each data type analysis
    for j=1:size(datasets,1)
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
            % get the features from each analysis
            for k=1:length(folders)
                data=readInputTableForPipeline([folders{k} filesep 'feature_importance' filesep 'final_feature_importance.csv']);
                table{1,cnt} = [datasets{j,1} ' ' lower(datatypes{i,1}) ' ' lower(folders{k})];
                table(2:size(data,1),cnt) = data(2:end,1);
                table(2:size(data,1),cnt+1) = data(2:end,2);
                cnt=cnt+2;
            end
        end
        cd(rootDir)
    end
end
writetable(cell2table(table),[rootDir filesep 'results' filesep 'strains' filesep 'Table_S6.csv'],'WriteVariableNames',false);

