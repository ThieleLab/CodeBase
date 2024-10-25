% get the distribution of reactions by subsystem across phyla in all
% reconstructions

loadedData=readInputTableForPipeline(['Model_properties_analysis' filesep '150k_90k_combined_text_files' filesep 'ReactionPresence_150k_90k_combined_refined.txt']);

phyla=cell(size(loadedData,1),1);

% get phylum information for each microbe
info=readInputTableForPipeline(['Model_properties_analysis' filesep '150k' filesep 'SequencedGenomesInfo.txt']);
phylumCol=find(strcmp(info(1,:),'Phylum'));
for i=2:size(loadedData,1)
    findStrain=find(strcmp(info(:,1),loadedData{i,1}));
    if ~isempty(findStrain)
        phyla{i,1}=info{findStrain,phylumCol};
    end
    if isempty(phyla{i,1})
        phyla{i,1}='N/A';
    end
end

info=readInputTableForPipeline(['Model_properties_analysis' filesep '90k' filesep 'mags-gut_qs50_checkm_adaptedToPipeline.txt']);
phylumCol=find(strcmp(info(1,:),'Phylum'));
for i=2:size(loadedData,1)
    findStrain=find(strcmp(info(:,1),loadedData{i,1}));
    if ~isempty(findStrain)
        phyla{i,1}=info{findStrain,phylumCol};
    end
    if isempty(phyla{i,1})
        phyla{i,1}='N/A';
    end
end

% get all unique phyla
uniquePhyla=unique(phyla(2:end));

% now get the number of phylum representatives for each reaction
data=zeros(length(uniquePhyla),size(loadedData,2)-1);
for i=1:length(uniquePhyla)
    getPhylum=find(strcmp(phyla(:,1),uniquePhyla{i}));
    for k=1:length(getPhylum)
        for j=2:size(loadedData,2)
            % fix some temporary empty spaces in the data
            if isempty(loadedData{getPhylum(k),j})
                loadedData{getPhylum(k),j}='0';
            end
            data(i,j-1)=data(i,j-1)+str2double(loadedData{getPhylum(k),j});
        end
    end
end
save(['Model_properties_analysis' filesep 'Summary_for_figures' filesep 'ReactionsByPhylum.mat'],'data');

% get the subsystems
database=loadVMHDatabase;

subs=cell(size(loadedData,2)-1,1);
for i=2:size(loadedData,2)
    findRxn=find(strcmp(database.reactions(:,1),loadedData{1,i}));
    if ~isempty(findRxn)
        subs{i-1,1}=database.reactions{findRxn,12};
    end
end

% summarize reaction distribution by phylum by subsystem
uniqueSubs=unique(subs);
dataSub=zeros(length(uniquePhyla),length(uniqueSubs));
for i=1:size(data,2)
    getSub=find(strcmp(uniqueSubs,subs{i,1}));
    for j=1:size(data,1)
        dataSub(j,getSub)=data(j,getSub)+data(j,i);
    end
end

% export as table
table=cell(length(uniquePhyla)+1,length(uniqueSubs)+1);
table(2:end,1)=uniquePhyla;
table(1,2:end)=uniqueSubs;
table(2:end,2:end)=num2cell(dataSub);
writetable(cell2table(table),['Model_properties_analysis' filesep 'Summary_for_figures' filesep 'PhylaBySubsystem.txt'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
