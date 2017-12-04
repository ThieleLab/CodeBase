% This script creates microbe-microbe models in every combination for all 
% 773  AGORA models. These pairwise modelsa were used in Magnusdottir et
% al., Nat Biotech 2017, to simulate pairwise interactions (results shown
% in Figure 5).

% NOTE: this process is very time-consuming due to the large number of 
% models created. If you are only interested in certain pairs, creating a
% reduced input file containing only the desired microbe models is
% recommended.

% To simulate the pairwise interactions afterwards, please see the
% functions already implemented in Cobra Toolbox v.3.
% cobratoolbox\src\analysis\multiSpecies
% cobratoolbox\papers\2016_pairwiseModeling_HumanGutMicrobiota
% Execute the script simulationPairwiseInteractions with the list of
% pairwise models created in this script (file "pairedModelsList") as the input.

% Almut Heinken, 01.12.2017

% initialize the COBRA Toolbox
initCobraToolbox

% Import file with information on AGORA including reconstruction names
currentDir = pwd;
[~,infoFile,~]=xlsread(strcat(currentDir,'\InputFiles\','ModelInformation.xlsx'));

% create a new folder where the pairwise models are saved
mkdir(currentDir,'pairwiseModels')

pairedModelsList{1, 1} = 'pairedModelID';
pairedModelsList{1, 2} = 'ModelID1';
pairedModelsList{1, 3} = 'Strain1';
pairedModelsList{1, 4} = 'Biomass1';
pairedModelsList{1, 5} = 'ModelID2';
pairedModelsList{1, 6} = 'Strain2';
pairedModelsList{1, 7} = 'Biomass2';
modelList=2;
for i=2:size(infoFile,1)
    load(strcat(currentDir,'\AGORA\',infoFile{i,3},'.mat'));
    model1=model;
    bioID1=find(strncmp(model1.rxns,'biomass',7));
    % make sure pairs are only generated once
    % and models aren't paired with themselves
    for j=i+1:length(infoFile)
        load(strcat(currentDir,'\AGORA\',infoFile{j,3},'.mat'));
        model2=model;
        bioID2=find(strncmp(model2.rxns,'biomass',7));
        % script "createMultiSpeciesModel" will join the two microbes
        % need to create file with the two models and two corresponding
        % name tags as input for createMultipleSpeciesModel
        models{1, 1} = model1;
        models{2, 1} = model2;
        nameTagsModels{1, 1} = strcat(infoFile{i,3}, '_');
        nameTagsModels{2, 1} = strcat(infoFile{j,3}, '_');
        % use the function createMultipleSpeciesModel with the required
        % inputs
        [pairedModel] = createMultipleSpeciesModel(models, nameTagsModels);

        % create coupling constraints: the reactions of each individual
        % microbe need to be coupled to its own biomass objective function.
        % For this, it is necessary to find all reactions in each paired
        % microbe by looking for the name tags each microbe received (see
        % "nameTagsModels" above). The biomass objective function is also
        % found by retrieving the name tag + name of the biomass objective
        % function. The coupling factor c here is 400 with a threshold u of
        % 0, this may be edited as convenient.
        
        [pairedModel]=coupleRxnList2Rxn(pairedModel,pairedModel.rxns(strmatch(strcat(infoFile{i,3},'_'),pairedModel.rxns)),strcat(infoFile{i,3},'_',model1.rxns(bioID1)),400,0);
        [pairedModel]=coupleRxnList2Rxn(pairedModel,pairedModel.rxns(strmatch(strcat(infoFile{j,3},'_'),pairedModel.rxns)),strcat(infoFile{j,3},'_',model2.rxns(bioID2)),400,0);
        
        save(strcat(currentDir,'\pairwiseModels\pairedModel_',infoFile{i,3},'_',infoFile{j,3},'.mat'),'pairedModel');
        % keep track of the generated models and populate the output file with
        % pairwise model information
        pairedModelsList{modelList,1}=strcat('pairedModel','_',infoFile{i,3},'_',infoFile{j,3},'.mat');
        pairedModelsList{modelList,2}=infoFile{i,3};
        pairedModelsList{modelList,3}=infoFile{i,1};
        pairedModelsList{modelList,4}=model1.rxns(bioID1);
        pairedModelsList{modelList,5}=infoFile{j,3};
        pairedModelsList{modelList,6}=infoFile{j,1};
        pairedModelsList{modelList,7}=model2.rxns(bioID2);
        modelList=modelList+1;
        if j==774
            save(strcat(currentDir,'\pairwiseModels\','pairedModelsList'),'pairedModelsList');
        end
    end
end
save(strcat(currentDir,'\pairwiseModels\','pairedModelsList'),'pairedModelsList');
