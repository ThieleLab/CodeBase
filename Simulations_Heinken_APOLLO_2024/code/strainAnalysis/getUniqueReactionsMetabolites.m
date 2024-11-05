
% get unique reactions and metabolites for Venn diagrams (Figures 2f, S1)

clear all
rootDir = pwd;

uniqueRxns = {
    'Almeida','Pasolli','APOLLO','AGORA2'
};

uniqueMets = {
    'Almeida','Pasolli','APOLLO','AGORA2'
};

allRxns = {};
allMets = {};

% Almeida
loadedData = parquetread([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'Pasolli_Almeida_parquet_files' filesep 'reactionPresence_Almeida_refined.parquet']);
reactions = loadedData.Properties.VariableNames;
reactions = regexprep(reactions,'x','','once');
metabolites = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Almeida' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_Almeida_refined.txt']);

uniqueRxns(2:length(reactions)+1,1) = reactions;
uniqueMets(2:size(metabolites,2),1) = metabolites(1,2:end);
allRxns = union(allRxns,reactions);
allMets = union(allMets,metabolites(1,2:end));

% Pasolli
loadedData = parquetread([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'Pasolli_Almeida_parquet_files' filesep 'reactionPresence_Pasolli_refined.parquet']);
reactions = loadedData.Properties.VariableNames;
reactions = regexprep(reactions,'x','','once');
metabolites = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Pasolli' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_Pasolli_refined.txt']);

uniqueRxns(2:length(reactions)+1,2) = reactions;
uniqueMets(2:size(metabolites,2),2) = metabolites(1,2:end);
allRxns = union(allRxns,reactions);
allMets = union(allMets,metabolites(1,2:end));

uniqueRxns(2:size(allRxns,1)+1,3) = allRxns;
uniqueMets(2:size(allMets,1)+1,3) = allMets;

% AGORA2
reactions = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_AGORA2' filesep 'ReactionMetabolitePresence' filesep 'ReactionPresence_AGORA2_refined.txt']);
metabolites = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_AGORA2' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_AGORA2_refined.txt']);

uniqueRxns(2:size(reactions,2),4) = reactions(1,2:end)';
uniqueMets(2:size(metabolites,2),4) = metabolites(1,2:end)';

writetable(cell2table(uniqueRxns),[rootDir filesep 'results' filesep 'strains' filesep 'UniqueReactions.csv'],'writeVariableNames',false)
writetable(cell2table(uniqueMets),[rootDir filesep 'results' filesep 'strains' filesep 'UniqueMetabolites.csv'],'writeVariableNames',false)

