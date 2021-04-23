% This script calculates absolute reaction presence, reaction abundance,
% and subsystem abundance in the personalized microbiome models generated 
% for the PLEASE/COMBO cohort.
% Almut Heinken, 03/21

%% Get the reactions to compute
abundance = readtable(abunFilePath, 'ReadVariableNames', false);
abundance = table2cell(abundance);
if isnumeric(abundance{2, 1})
    abundance(:, 1) = [];
end

% Do not include biomass and exchange/demand-not very informative
reactions2Compute = {};
    for i = 2:size(abundance, 1)
        load([agoraPath filesep abundance{i, 1} '.mat']);
        %     model = readCbModel([modelPath filesep abundance{i, 1} '.mat']);
        modelsList{i, 1} = model;
    end
    % get model list from abundance input file
    for i = 2:size(abundance, 1)
        model = modelsList{i, 1};
        reactions2Compute = vertcat(model.rxns, reactions2Compute);
    end
    reactions2Compute = unique(reactions2Compute);
    % Remove biomass and exchange/demand reactions
    delRows=[];
    cnt=1;
    for i=2:size(reactions2Compute,1)
        if strncmp(reactions2Compute{i,1},'biomass',7) || strncmp(reactions2Compute{i,1},'EX_',3) || strncmp(reactions2Compute{i,1},'DM_',3) || strncmp(reactions2Compute{i,1},'sink_',5)
            delRows(cnt,1)=i;
            cnt=cnt+1;
        end
    end
    reactions2Compute(delRows,:)=[];

%% Calculate reaction presence/absence
[ReactionPresence,ReactionPresenceDifferent]=calculateReactionPresence(abunFilePath, agoraPath,reactions2Compute);
cell2csv([rxnsPath filesep 'ReactionPresence.csv'],ReactionPresenceDifferent);

%% Calculate reaction abundance
ReactionAbundance=calculateReactionAbundance(abunFilePath,agoraPath,'AGORA_infoFile.xlsx',reactions2Compute,numWorkers);

% export the files as csv sheets
cell2csv([rxnsPath filesep 'ReactionAbundance.csv'],ReactionAbundance.Total');
cell2csv([rxnsPath filesep 'ReactionAbundance_Phylum.csv'],ReactionAbundance.Phylum');
cell2csv([rxnsPath filesep 'ReactionAbundance_Genus.csv'],ReactionAbundance.Genus');

%% Calculate subsystem abundance

reactionAbundancePath = [rxnsPath filesep 'ReactionAbundance.csv'];
SubsystemAbundance = calculateSubsystemAbundance(reactionAbundancePath);
SubsystemAbundance(:,1)=strrep(SubsystemAbundance(:,1),',',' ');
cell2csv([rxnsPath filesep 'SubsystemAbundance.csv'],SubsystemAbundance);
