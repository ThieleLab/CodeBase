% get reconstruction stats and plot them separated by taxon, body site,
% continent, lifestyle, and age group

load(['Manuscript' filesep 'Computation_Figure_2' filesep 'data_combined.mat'])

% get the names for each reconstruction
% Almeida reconstructions
path_90kModels='/Users/almut.heinken/OneDrive - National University of Ireland, Galway/90000_genomes_Almeida_2019/curatedReconstructions';
dInfo = dir(path_90kModels);
modelList={dInfo.name};
modelList=modelList';
modelList(find(~contains(modelList,'.mat')),:)=[];
models=strrep(modelList,'.mat','');

% Pasolli reconstructions
path_150kModels='/Users/almut.heinken/OneDrive - National University of Ireland, Galway/150k_Project_Models/curatedReconstructions';
dInfo = dir(path_150kModels);
modelList={dInfo.name};
modelList=modelList';
modelList(find(~contains(modelList,'.mat')),:)=[];
models=vertcat(models,strrep(modelList,'.mat',''));

% read in strain information
info=readInputTableForPipeline(['Model_properties_analysis' filesep '90k_150k_combined' filesep 'Combined_taxonomy_info.xlsx']);
for i=2:size(info,1)
    for j=2:size(info,2)
        if ~ischar(info{i,j}) || isempty(info{i,j})
        info{i,j}='';
        end
    end
end
% cd(['Model_properties_analysis' filesep 'Summary_for_figures'])
% mkdir('Strain_properties_by_group')
% cd('Strain_properties_by_group')

datatypes={'Number of reactions','Number of metabolites','Number of genes','Aerobic growth','Anaerobic growth','Aerobic ATP production','Anaerobic ATP production'};

% export the data for top 30 entries in each category (if taxa) or
% otherwise all groups for each type of data (number of reactions, number
% of metabolites...)
for i=2:size(info,2)
    i
    [feats,~,J]=unique(info(2:end,i));
    if length(feats)>30
        cutoff=30;
    else
        cutoff=length(feats);
    end
    cnt = histc(J, 1:numel(feats));

    % removed unspecified taxa
    TF = endsWith(feats,' sp');
    cnt(find(TF==true),:)=[];
    feats(find(TF==true),:)=[];
    TF = endsWith(feats,' unclassified');
    cnt(find(TF==true),:)=[];
    feats(find(TF==true),:)=[];
    TF = endsWith(feats,' bacterium');
    cnt(find(TF==true),:)=[];
    feats(find(TF==true),:)=[];
    cnt(find(strcmp(feats,'')),:)=[];
    feats(find(strcmp(feats,'')),:)=[];
    cnt(find(strcmp(feats,'N/A')),:)=[];
    feats(find(strcmp(feats,'N/A')),:)=[];
    cnt(find(strcmp(feats,'NA')),:)=[];
    feats(find(strcmp(feats,'NA')),:)=[];
    
    [A,I]=sort(cnt,'descend');
    labels=feats(I(1:cutoff),1);

    % print the data
    for j=1:length(datatypes)
        data = {};
        for k=1:length(labels)
            data{1,k} = labels{k};
            findInd=find(strcmp(info(:,i),labels{k}));
            for l=1:length(findInd)
                findModel = find(strcmp(models,info{findInd(l),1}));
                data{l+1,k} = data_combined(findModel,j);
            end
        end
        cell2csv(['Model_properties_analysis' filesep 'Summary_for_figures' filesep 'Strain_properties_by_group' filesep info{1,i} '_' strrep(datatypes{j},' ','_') '.csv'],data)
    end
end
