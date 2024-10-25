% plot random assembly of reconstructions and how many of total reactions
% are added with each new reconstruction

%% Reaction presence
% all reconstructions
loadedData=readInputTableForPipeline(['Model_properties_analysis' filesep '90k_properties' filesep 'ReactionMetabolitePresence' filesep 'ReactionPresence_90k_refined.txt']);
data=loadedData(2:end,2:end);
% get a random order of reconstructions
p = randperm(size(data,1));

totalRxns=[];
indModels=[];
sumRxns=[];
for i=1:length(p)
    presentRxns=find(strcmp(data(p(i),:),'1'));
    totalRxns=union(totalRxns,presentRxns);
    indModels(i,1)=i;
    sumRxns(i,1)=length(totalRxns);
end

% plot the number of reactions added with each models
figure
scatter(indModels,sumRxns,'.','k')
title('Reaction distribution in 90k reconstructions')
xlabel('Number of reconstructions added')
ylabel('Number of reactions added')
print(['Analysis_microbiome_models' filesep 'Analysis_Plots' filesep 'Reaction_distribution_90k'],'-dpng','-r300')

% 150k
loadedData=readInputTableForPipeline(['Model_properties_analysis' filesep '150k_properties' filesep 'ReactionMetabolitePresence' filesep 'ReactionPresence_150k_refined.txt']);
data=loadedData(2:end,2:end);
% get a random order of reconstructions
p = randperm(size(data,1));

totalRxns=[];
indModels=[];
sumRxns=[];
for i=1:length(p)
    presentRxns=find(strcmp(data(p(i),:),'1'));
    totalRxns=union(totalRxns,presentRxns);
    indModels(i,1)=i;
    sumRxns(i,1)=length(totalRxns);
end

% plot the number of reactions added with each models
figure
scatter(indModels,sumRxns,'.','k')
title('Reaction distribution in 150k reconstructions')
xlabel('Number of reconstructions added')
ylabel('Number of reactions added')
print(['Analysis_microbiome_models' filesep 'Analysis_Plots' filesep 'Reaction_distribution_150k'],'-dpng','-r300')

% both datasets combined
loadedData=readInputTableForPipeline(['Model_properties_analysis' filesep '90k_150k_combined' filesep '150k_90k_combined_text_files' filesep 'ReactionPresence_150k_90k_combined_refined.txt']);
data=loadedData(2:end,2:end);
% get a random order of reconstructions
p = randperm(size(data,1));

totalRxns=[];
indModels=[];
sumRxns=[];
for i=1:length(p)
    presentRxns=find(strcmp(data(p(i),:),'1'));
    totalRxns=union(totalRxns,presentRxns);
    indModels(i,1)=i;
    sumRxns(i,1)=length(totalRxns);
end
% get all cases where each reaction has been added
% allRxns=find(sumRxns==max(sumRxns));
% indModels(allRxns(2:end),:)=[];
% sumRxns(allRxns(2:end),:)=[];

% plot the number of reactions added with each models
figure
scatter(indModels,sumRxns,'.','k')
title('Reaction distribution in all reconstructions')
xlabel('Number of reconstructions added')
ylabel('Number of reactions added')
print(['Analysis_microbiome_models' filesep 'Analysis_Plots' filesep 'Reaction_distribution_combined'],'-dpng','-r300')

%% Metabolite presence
% 90k
loadedData=readInputTableForPipeline(['Model_properties_analysis' filesep '90k_properties' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_90k_refined.txt']);
data=loadedData(2:end,2:end);
% get a random order of reconstructions
p = randperm(size(data,1));

totalRxns=[];
indModels=[];
sumRxns=[];
for i=1:length(p)
    presentRxns=find(strcmp(data(p(i),:),'1'));
    totalRxns=union(totalRxns,presentRxns);
    indModels(i,1)=i;
    sumRxns(i,1)=length(totalRxns);
end

% plot the number of reactions added with each models
figure
scatter(indModels,sumRxns,'.','k')
title('Metabolite distribution in 90k reconstructions')
xlabel('Number of reconstructions added')
ylabel('Number of reactions added')
print(['Analysis_microbiome_models' filesep 'Analysis_Plots' filesep 'Metabolite_distribution_90k'],'-dpng','-r300')

% 150k
loadedData=readInputTableForPipeline(['Model_properties_analysis' filesep '150k_properties' filesep 'ReactionMetabolitePresence' filesep 'MetabolitePresence_150k_refined.txt']);
data=loadedData(2:end,2:end);
% get a random order of reconstructions
p = randperm(size(data,1));

totalRxns=[];
indModels=[];
sumRxns=[];
for i=1:length(p)
    presentRxns=find(strcmp(data(p(i),:),'1'));
    totalRxns=union(totalRxns,presentRxns);
    indModels(i,1)=i;
    sumRxns(i,1)=length(totalRxns);
end

% plot the number of reactions added with each models
figure
scatter(indModels,sumRxns,'.','k')
title('Metabolite distribution in 150k reconstructions')
xlabel('Number of reconstructions added')
ylabel('Number of reactions added')
print(['Analysis_microbiome_models' filesep 'Analysis_Plots' filesep 'Metabolite_distribution_150k'],'-dpng','-r300')

% both datasets combined
loadedData=readInputTableForPipeline(['Model_properties_analysis' filesep '90k_150k_combined' filesep '150k_90k_combined_text_files' filesep 'MetabolitePresence_150k_90k_combined_refined.txt']);
data=loadedData(2:end,2:end);
% get a random order of reconstructions
p = randperm(size(data,1));

totalRxns=[];
indModels=[];
sumRxns=[];
for i=1:length(p)
    presentRxns=find(strcmp(data(p(i),:),'1'));
    totalRxns=union(totalRxns,presentRxns);
    indModels(i,1)=i;
    sumRxns(i,1)=length(totalRxns);
end

% plot the number of reactions added with each models
figure
scatter(indModels,sumRxns,'.','k')
title('Metabolite distribution in all reconstructions')
xlabel('Number of reconstructions added')
ylabel('Number of reactions added')
print(['Analysis_microbiome_models' filesep 'Analysis_Plots' filesep 'Metabolite_distribution_combined'],'-dpng','-r300')
