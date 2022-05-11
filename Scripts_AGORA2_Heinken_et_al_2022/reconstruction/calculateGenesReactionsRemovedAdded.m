                                                                                                                                                        
% Get the number of genes that were added/removed for each draft
% reconstruction during refinement.

dInfo = dir(refinedFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

geneCount={'Model_ID','Genes_only_in_draft','Genes_only_in_curated','Genes_in_both'};
rxnCount={'Model_ID','Reactions_only_in_draft','Reactions_only_in_curated','Reactions_in_both'};

% Load reaction and metabolite database
database=loadVMHDatabase;

for i=1:length(modelList)
    i
    geneCount{i+1,1}=strrep(modelList{i},'.mat','');
    rxnCount{i+1,1}=strrep(modelList{i},'.mat','');
    draftModel=readCbModel([translatedDraftsFolder filesep modelList{i}]);
    % translate draft model to get the right numbers
    biomassReaction=draftModel.rxns{find(strncmp(draftModel.rxns,'bio',3)),1};
    draftModel = translateKBaseModel2VMHModel(draftModel, biomassReaction,database);
    
    curatedModel=readCbModel([refinedFolder filesep modelList{i}]);
    % genes
    [C,IA]=setdiff(draftModel.genes,curatedModel.genes);
    geneCount{i+1,2}=length(C);
    [C,IA]=setdiff(curatedModel.genes,draftModel.genes);
    geneCount{i+1,3}=length(C);
    [C,IA]=intersect(draftModel.genes,curatedModel.genes);
    geneCount{i+1,4}=length(C);
    % reactions
    [C,IA]=setdiff(draftModel.rxns,curatedModel.rxns);
    rxnCount{i+1,2}=length(C);
    [C,IA]=setdiff(curatedModel.rxns,draftModel.rxns);
    rxnCount{i+1,3}=length(C);
    [C,IA]=intersect(draftModel.rxns,curatedModel.rxns);
    rxnCount{i+1,4}=length(C);
end

writetable(cell2table(geneCount),[propertiesFolder filesep 'Genes_Removed_Added.xlsx'],'FileType','spreadsheet','WriteVariableNames',false);
writetable(cell2table(rxnCount),[propertiesFolder filesep 'Reactions_Removed_Added.xlsx'],'FileType','spreadsheet','WriteVariableNames',false);
