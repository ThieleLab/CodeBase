
% get unique reactions and metabolites for Venn diagrams

uniqueRxns = {
    'Almeida','Pasolli','APOLLO','AGORA2'
};

uniqueMets = {
    'Almeida','Pasolli','APOLLO','AGORA2'
};

allRxns = {};
allMets = {};

% Almeida
reactions = readInputTableForPipeline(['Model_properties_analysis' filesep '90k_properties' filesep 'ReactionMetabolitePresence' filesep 'ReactionPresence_90k_refined.txt']);
metabolites = readInputTableForPipeline(['Model_properties_analysis' filesep '90k_properties' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_90k_refined.txt']);

uniqueRxns(2:size(reactions,1),1) = reactions(2:end,1);
uniqueMets(2:size(metabolites,1),1) = metabolites(2:end,1);
allRxns = union(allRxns,reactions(2:end,1));
allMets = union(allMets,metabolites(2:end,1));

% Pasolli
reactions = readInputTableForPipeline(['Model_properties_analysis' filesep '150k_properties' filesep 'ReactionMetabolitePresence' filesep 'ReactionPresence_150k_refined.txt']);
metabolites = readInputTableForPipeline(['Model_properties_analysis' filesep '150k_properties' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_150k_refined.txt']);

uniqueRxns(2:size(reactions,1),2) = reactions(2:end,1);
uniqueMets(2:size(metabolites,1),2) = metabolites(2:end,1);
allRxns = union(allRxns,reactions(2:end,1));
allMets = union(allMets,metabolites(2:end,1));

uniqueRxns(2:size(allRxns,1)+1,3) = allRxns;
uniqueMets(2:size(allMets,1)+1,3) = allMets;

% AGORA2
reactions = readInputTableForPipeline(['Model_properties_analysis' filesep 'AGORA2' filesep 'ReactionMetabolitePresence' filesep 'ReactionPresence_AGORA2_refined.txt']);
metabolites = readInputTableForPipeline(['Model_properties_analysis' filesep 'AGORA2' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_AGORA2_refined.txt']);

uniqueRxns(2:size(reactions,1),4) = reactions(2:end,1);
uniqueMets(2:size(metabolites,1),4) = metabolites(2:end,1);

cell2csv('UniqueReactions.csv','uniqueRxns')
cell2csv('UniqueMetabolites.csv','uniqueMets')
