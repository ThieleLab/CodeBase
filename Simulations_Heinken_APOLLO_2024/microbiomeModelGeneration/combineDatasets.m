
%% first combine the datasets
filePaths={
    'Nasal_cavity','Results_nasal_cavity';
    'Vagina','Results_vagina';
    'Skin','Results_skin';
    'Gut','Results_gut';
    };

% combine files for all datasets

% combine model statistics

stats={'ModelID','Reactions','Metabolites'};
cnt=2;
for i=1:length(filePaths)
    data = table2cell(readtable([filePaths{i,2} filesep 'ModelStatistics.csv'], 'ReadVariableNames', false));
    for j=2:size(data,1)
        stats{cnt,1}=data{j,1};
        stats{cnt,2}=data{j,2};
        stats{cnt,3}=data{j,3};
        cnt=cnt+1;
    end
end
stats(:,1)=strrep(stats(:,1),',',' ');
cell2csv('ModelStatisticsCombined.csv',stats)

% combine reaction abundances
reactions={''};
for i=1:length(filePaths)
    data = readInputTableForPipeline([filePaths{i,2} filesep 'ReactionAbundance.csv']);
    reactions=union(reactions,data(2:end,1));
end

cnt=2;
for i=1:length(filePaths)
    data = readInputTableForPipeline([filePaths{i,2} filesep 'ReactionAbundance.csv']);
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
cell2csv('ReactionAbundanceCombined.csv',reactions)

% calculate the combined subsystem abundance

subsystemAbundance = calculateSubsystemAbundance('ReactionAbundanceCombined.csv');
subsystemAbundance(:,1)=strrep(subsystemAbundance(:,1),',',' ');
cell2csv('SubsystemAbundanceCombined.csv',subsystemAbundance')

cell2csv('ReactionAbundanceCombined.csv',reactions')

% calculate reaction presence
reactions=readInputTableForPipeline('ReactionAbundanceCombined.csv');
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
cell2csv('ReactionPresenceCombined.csv',reactions)

% combine organism abundances

filePaths={
    'normalized_gut_abundances.csv'
    'normalized_nasal_cavity_abundances.csv';
    'normalized_skin_abundances.csv'
    'normalized_vagina_abundances.csv';
    };

abundance={''};
for i=1:length(filePaths)
    data = readInputTableForPipeline(['abundances_from_metagenomes' filesep filePaths{i,1}]);
    abundance=union(abundance,data(2:end,1));
end

cnt=2;
for i=1:length(filePaths)
    data = readInputTableForPipeline(['abundances_from_metagenomes' filesep filePaths{i,1}]);
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
cell2csv('normalizedAbundanceCombined.csv',abundance)
