
% get unique reactions and metabolites for Venn diagrams (Figures 2f, S1)

uniqueRxns = {
    'Almeida','Pasolli','APOLLO','AGORA2'
};

uniqueMets = {
    'Almeida','Pasolli','APOLLO','AGORA2'
};

allRxns = {};
allMets = {};

% Almeida
reactions = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Almeida' filesep 'ReactionMetabolitePresence' filesep 'ReactionPresence_Almeida_refined.txt']);
metabolites = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Almeida' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_Almeida_refined.txt']);

uniqueRxns(2:size(reactions,1),1) = reactions(2:end,1);
uniqueMets(2:size(metabolites,1),1) = metabolites(2:end,1);
allRxns = union(allRxns,reactions(2:end,1));
allMets = union(allMets,metabolites(2:end,1));

% Pasolli
reactions = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Pasolli' filesep 'ReactionMetabolitePresence' filesep 'ReactionPresence_Pasolli_refined.txt']);
metabolites = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Pasolli' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_Pasolli_refined.txt']);

uniqueRxns(2:size(reactions,1),2) = reactions(2:end,1);
uniqueMets(2:size(metabolites,1),2) = metabolites(2:end,1);
allRxns = union(allRxns,reactions(2:end,1));
allMets = union(allMets,metabolites(2:end,1));

uniqueRxns(2:size(allRxns,1)+1,3) = allRxns;
uniqueMets(2:size(allMets,1)+1,3) = allMets;

% AGORA2
reactions = readInputTableForPipeline([rootDir filesep 'data' filesep 'refinedModelProperties_AGORA2' filesep 'ReactionMetabolitePresence' filesep 'ReactionPresence_AGORA2_refined.txt']);
metabolites = readInputTableForPipeline([rootDir filesep 'data' filesep 'refinedModelProperties_AGORA2' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_AGORA2_refined.txt']);

uniqueRxns(2:size(reactions,1),4) = reactions(2:end,1);
uniqueMets(2:size(metabolites,1),4) = metabolites(2:end,1);

cell2csv([rootDir filesep 'results' filesep 'strains' filesep 'UniqueReactions.csv'],'uniqueRxns')
cell2csv([rootDir filesep 'results' filesep 'strains' filesep 'UniqueMetabolites.csv'],'uniqueMets')
