
%% get the gut samples from 90k dataset
metadata90kPath='C:\Users\Almut Heinken\OneDrive - National University of Ireland, Galway\90000_genomes_Almeida_2019\MicrobiomeModelling\gut_metagenomic_samples_metadata.csv';
rxnAbun90kPath='C:\Users\Almut Heinken\OneDrive - National University of Ireland, Galway\90000_genomes_Almeida_2019\MicrobiomeModelling\Reaction_abundances\ReactionAbundance.csv';
subAbun90kPath='C:\Users\Almut Heinken\OneDrive - National University of Ireland, Galway\90000_genomes_Almeida_2019\MicrobiomeModelling\Reaction_abundances\SubsystemAbundance.csv';

%% first combine the datasets
filePaths={
    'Nasal_cavity';
    'Vagina';
    'Skin';
    };

% combine files for all datasets

% create a metadata file
metadata={'Sample','Body site'};
cnt=2;
for i=1:length(filePaths)
    data = table2cell(readtable([filePaths{i,1} filesep 'SubsystemAbundance.csv'], 'ReadVariableNames', false));
    for j=2:size(data,2)
        metadata{cnt,1}=data{1,j};
        metadata{cnt,2}=filePaths{i};
        cnt=cnt+1;
    end
end

% add the gut samples from 90k
data = table2cell(readtable(metadata90kPath, 'ReadVariableNames', false));
% only US to be comparable
IA=find(~strcmp(data(:,14),'United States'));
data(IA(2:end),:)=[];
% remove infant and elderly samples
IA=find(ismember(data(:,7),{'Infant','Elderly'}));
data(IA(2:end),:)=[];
% remove disease cases
IA=find(~ismember(data(:,5),{'Healthy','NA'}));
data(IA(2:end),:)=[];
for j=2:size(data,1)
    metadata{cnt,1}=data{j,1};
    metadata{cnt,2}='Gut';
    cnt=cnt+1;
end
cell2csv('metadata.csv',metadata)

% combine reaction abundances
reactions={''};
for i=1:length(filePaths)
    data = table2cell(readtable([filePaths{i,1} filesep 'ReactionAbundance.csv'], 'ReadVariableNames', false));
    reactions=union(reactions,data(2:end,1));
end
data = table2cell(readtable(rxnAbun90kPath, 'ReadVariableNames', false));
reactions=union(reactions,data(2:end,1));

cnt=2;
for i=1:length(filePaths)
    data = table2cell(readtable([filePaths{i,1} filesep 'ReactionAbundance.csv'], 'ReadVariableNames', false));
    for j=2:size(data,2)
        reactions{1,cnt}=data{1,j};
        reactions(2:end,cnt)={'0'};
        for k=2:size(data,1)
            rxnInd=find(strcmp(reactions(:,1),data{k,1}));
            reactions{rxnInd,cnt}=data{k,j};
        end
        cnt=cnt+1;
    end
end

% add the gut samples
data = table2cell(readtable(rxnAbun90kPath, 'ReadVariableNames', false));
[C,IA]=setdiff(data(1,:),metadata(:,1),'stable');
data(:,IA(2:end))=[];
for j=2:size(data,2)
    reactions{1,cnt}=data{1,j};
    reactions(2:end,cnt)={'0'};
    for k=2:size(data,1)
        rxnInd=find(strcmp(reactions(:,1),data{k,1}));
        reactions{rxnInd,cnt}=data{k,j};
    end
    cnt=cnt+1;
end
cell2csv('ReactionAbundanceCombined.csv',reactions)

% combine subsystem abundances
reactions={''};
for i=1:length(filePaths)
    data = table2cell(readtable([filePaths{i,1} filesep 'SubsystemAbundance.csv'], 'ReadVariableNames', false));
    reactions=union(reactions,data(2:end,1));
end
data = table2cell(readtable(subAbun90kPath, 'ReadVariableNames', false));
reactions=union(reactions,data(1,2:end));

cnt=2;
for i=1:length(filePaths)
    data = table2cell(readtable([filePaths{i,1} filesep 'SubsystemAbundance.csv'], 'ReadVariableNames', false));
    for j=2:size(data,2)
        reactions{1,cnt}=data{1,j};
        reactions(2:end,cnt)={'0'};
        for k=2:size(data,1)
            rxnInd=find(strcmp(reactions(:,1),data{k,1}));
            reactions{rxnInd,cnt}=data{k,j};
        end
        cnt=cnt+1;
    end
end

% add the gut samples
data = table2cell(readtable(subAbun90kPath, 'ReadVariableNames', false));
data=data';
[C,IA]=setdiff(data(1,:),metadata(:,1),'stable');
data(:,IA(2:end))=[];
for j=2:size(data,2)
    reactions{1,cnt}=data{1,j};
    reactions(2:end,cnt)={'0'};
    for k=2:size(data,1)
        rxnInd=find(strcmp(reactions(:,1),data{k,1}));
        reactions{rxnInd,cnt}=data{k,j};
    end
    cnt=cnt+1;
end
reactions(:,1)=strrep(reactions(:,1),',',' ');
cell2csv('SubsystemAbundanceCombined.csv',reactions)

%% then analyze the combined datsets

combinedDatasets={
    [pwd filesep 'ReactionAbundanceCombined.csv'],'Reaction abundance'
    [pwd filesep 'SubsystemAbundanceCombined.csv'],'Subsystem abundance'
    };

mkdir('ViolinPlots')
cd('ViolinPlots')

for i=2%1:length(combinedDatasets)
    data = table2cell(readtable(combinedDatasets{i,1}, 'ReadVariableNames', false));
    makeViolinPlots(data, metadata, 'plottedFeature',combinedDatasets{i,2},'unit','Relative abundance')
end

cd ..

mkdir('StatisticalAnalysis')
cd('StatisticalAnalysis')

for i=1:length(combinedDatasets)
    data = table2cell(readtable(combinedDatasets{i,1}, 'ReadVariableNames', false));
    [Statistics,significantFeatures] = performStatisticalAnalysis(data', metadata);
    Statistics(:,2)=strrep(Statistics(:,2),',',' ');
    cell2csv(['Statistics_' strrep(combinedDatasets{i,2}, ' ','_') '.csv'],Statistics)
    cell2csv(['SignificantFeatures_' strrep(combinedDatasets{i,2}, ' ','_') '.csv'],significantFeatures)
end

cd ..
