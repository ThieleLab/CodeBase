% This script creates microbe-microbe models in every combination for all 
% 773  AGORA models.
% Almut Heinken 2016-2017

% NOTE: this process is very time-consuming due to the large number of 
% models created. If you are only interested in certain pairs, creating a
% reduced input file containing only the desired microbe models is
% recommended.

% initialize the COBRA Toolbox
initCobraToolbox

% Import file with information on AGORA including reconstruction names
[~,labels,~]=xlsread('ModelInformation.xlsx');

currentDir = pwd;

pairedModelsList{1,1}='pairedModelID';
pairedModelsList{1,2}='ModelID1';
pairedModelsList{1,3}='ModelID2';
modelList=2;
for i=2:size(labels,1)
    load(strcat(currentDir,'\labels{i,3}','.mat'));
    % name the first model "model1" by default
    model1=model;
    % make sure pairs are only generated once
    % and models aren't paired with themselves
    for j=i+1:length(labels)
        load(strcat(currentDir,'\labels{j,3}','.mat'));
        % name the second model "model2" by default
        model2=model;
        modelsNameTags{1,1}=model1;
        modelsNameTags{2,1}=model2;
        modelsNameTags{1,2}=labels{i,3};
        modelsNameTags{2,2}=labels{j,3};
        % no host - host entry field empty
        [pairedModel] = createMultipleSpeciesModel([],[],2,modelsNameTags);
        
        [pairedModel]=coupleRxnList2Rxn(pairedModel,pairedModel.rxns(strmatch(labels{i,3},pairedModel.rxns)),strcat(labels{i,3},'biomass0'),400,0);
        [pairedModel]=coupleRxnList2Rxn(pairedModel,pairedModel.rxns(strmatch(labels{j,3},pairedModel.rxns)),strcat(labels{j,3},'biomass0'),400,0);
        
        save(strcat('Y:\Studies\Microbiome\Stefania\Microbiota_models\PairwiseModels_All_AGORA_Temp\pairedModel','_',labels{i,3},'_',labels{j,3},'.mat'),'pairedModel');
        % keep track of the generated models and populate the table with
        % information on created pairwise models
        pairedModelsList{modelList,1}=strcat('pairedModel','_',labels{i,3},'_',labels{j,3},'.mat');
        pairedModelsList{modelList,2}=labels{i,3};
        pairedModelsList{modelList,3}=labels{j,3};
        modelList=modelList+1;
        if j==774
            save pairedModelsList  pairedModelsList;
        end
    end
end
save pairedModelsList  pairedModelsList;
