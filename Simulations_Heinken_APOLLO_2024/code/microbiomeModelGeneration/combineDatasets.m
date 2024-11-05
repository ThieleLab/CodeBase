
% combine the microbiome model datasets from body sites and gut into one file

clear all
rootDir = pwd;

mkdir([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels'])

filePaths={
    ['BodySiteMicrobiomes' filesep 'Nasal_cavity']
    ['BodySiteMicrobiomes' filesep 'Vagina']
    ['BodySiteMicrobiomes' filesep  'Skin']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_1']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_2']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_3']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_4']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_5']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_6']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_7']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_8']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_9']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_10']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_11']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_12']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_13']
    ['GutMicrobiomes' filesep 'MicrobiomeModels_14']
    };

% combine model statistics

stats={'ModelID','Reactions','Metabolites'};
cnt=2;
for i=1:length(filePaths)
    data = table2cell(readtable([rootDir filesep 'data' filesep filePaths{i,1} filesep 'ModelStatistics.csv'], 'ReadVariableNames', false));
    for j=2:size(data,1)
        stats{cnt,1}=data{j,1};
        stats{cnt,2}=data{j,2};
        stats{cnt,3}=data{j,3};
        cnt=cnt+1;
    end
end
stats(:,1)=strrep(stats(:,1),',',' ');
writetable(cell2table(stats),[rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'ModelStatistics.csv'],'writeVariableNames',false)

% combine reaction abundances
reactions={''};
for i=1:length(filePaths)
    data = readInputTableForPipeline([filePaths{i,1} filesep 'ReactionAbundance.csv']);
    reactions=union(reactions,data(2:end,1));
end

cnt=2;
for i=1:length(filePaths)
    data = readInputTableForPipeline([filePaths{i,1} filesep 'ReactionAbundance.csv']);
    for j=2:size(data,2)
        j
        reactions{1,cnt}=data{1,j};
        reactions(2:end,cnt)={'0'};
        for k=2:size(data,1)
            rxnInd=find(strcmp(reactions(:,1),data{k,1}));
            reactions{rxnInd,cnt}=data{k,j};
        end
        cnt=cnt+1;
    end
end
reactions(:,1)=strrep(reactions(:,1),',',' ');
% remove biomass objective functions
bio=find(strncmp(reactions(:,1),'bio',3));
reactions(bio,:)=[];
writetable(cell2table(reactions),[rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'ReactionAbundance.csv'],'writeVariableNames',false)

% calculate the combined subsystem abundance
subsystemAbundance = calculateSubsystemAbundance([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'ReactionAbundance.csv']);
subsystemAbundance(:,1)=strrep(subsystemAbundance(:,1),',',' ');
writetable(cell2table(subsystemAbundance),[rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'SubsystemAbundance.csv'],'writeVariableNames',false)

% calculate reaction presence
reactions=readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'ReactionAbundance.csv']);
for i=2:size(reactions,1)
    for j=2:size(reactions,2)
        if ~isnumeric(reactions{i,j})
            reactions{i,j} = str2double(reactions{i,j});
        end
        if reactions{i,j} > 0.000001
            reactions{i,j} ='1';
        else
            reactions{i,j} ='0';
        end
    end
end
writetable(cell2table(reactions),[rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'ReactionPresence.csv'],'writeVariableNames',false)

% combine organism abundances

filePaths={
    'normalized_nasal_cavity_abundances.csv'
    'normalized_vagina_abundances.csv'
    'normalized_skin_abundances.csv'
    'normalized_gut_abundances.csv';
    };

abundance={''};
for i=1:length(filePaths)
    data = readInputTableForPipeline([rootDir filesep 'input' filesep filePaths{i,1}]);
    abundance=union(abundance,data(2:end,1));
end

cnt=2;
for i=1:length(filePaths)
    data = readInputTableForPipeline([rootDir filesep 'input' filesep filePaths{i,1}]);
    for j=2:size(data,2)
        j
        abundance{1,cnt}=data{1,j};
        abundance(2:end,cnt)={'0'};
        for k=2:size(data,1)
            rxnInd=find(strcmp(abundance(:,1),data{k,1}));
            abundance{rxnInd,cnt}=data{k,j};
        end
        cnt=cnt+1;
    end
end
writetable(cell2table(abundance),[rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'normalizedAbundance.csv'],'writeVariableNames',false)

