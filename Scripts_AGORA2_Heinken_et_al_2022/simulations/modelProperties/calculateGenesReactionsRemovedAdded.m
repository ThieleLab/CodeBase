
% Get the number of genes that were added/removed for each draft
% reconstruction during refinement.

resultsFolder = [rootDir filesep 'modelProperties' filesep];

dInfo = dir(refinedFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

geneCount={'Model_ID','Genes_only_in_draft','Genes_only_in_curated','Genes_in_both'};
rxnCount={'Model_ID','Reactions_only_in_draft','Reactions_only_in_curated','Reactions_in_both'};

% Load reaction and metabolite database
database=loadVMHDatabase;

% convert data types if necessary
database.metabolites(:,7)=strrep(database.metabolites(:,7),'NaN','');
database.metabolites(:,8)=strrep(database.metabolites(:,8),'NaN','');

for i=1:100:length(modelList)
    if length(modelList)-i>=steps-1
        endPnt=steps-1;
    else
        endPnt=length(modelList)-i;
    end

    geneCountTmp = {};
    rxnCountTmp = {};

    parfor j=i:i+endPnt
        draftModel=readCbModel([translatedDraftsFolder filesep modelList{j}]);
        % translate draft model to get the right numbers
        biomassReaction=draftModel.rxns{find(strncmp(draftModel.rxns,'bio',3)),1};
        draftModel = translateKBaseModel2VMHModel(draftModel, biomassReaction,database);

        curatedModel=readCbModel([refinedFolder filesep modelList{j}]);
        % genes
        [C,IA]=setdiff(draftModel.genes,curatedModel.genes);
        geneCountTmp{j}{1}=length(C);
        [C,IA]=setdiff(curatedModel.genes,draftModel.genes);
        geneCountTmp{j}{2}=length(C);
        [C,IA]=intersect(draftModel.genes,curatedModel.genes);
        geneCountTmp{j}{3}=length(C);
        % reactions
        [C,IA]=setdiff(draftModel.rxns,curatedModel.rxns);
        rxnCountTmp{j}{1}=length(C);
        [C,IA]=setdiff(curatedModel.rxns,draftModel.rxns);
        rxnCountTmp{j}{2}=length(C);
        [C,IA]=intersect(draftModel.rxns,curatedModel.rxns);
        rxnCountTmp{j}{3}=length(C);
    end

    for j=i:i+endPnt
        geneCount{j+1,1}=strrep(modelList{j},'.mat','');
        rxnCount{j+1,1}=strrep(modelList{j},'.mat','');
        for l=1:3
            geneCount{j+1,l+1}=geneCountTmp{j}{l};
            rxnCount{j+1,l+1}=rxnCountTmp{j}{l};
        end
    end
end

writetable(cell2table(geneCount),[resultsFolder filesep 'Genes_Removed_Added.xlsx'],'FileType','spreadsheet','WriteVariableNames',false);
writetable(cell2table(rxnCount),[resultsFolder filesep 'Reactions_Removed_Added.xlsx'],'FileType','spreadsheet','WriteVariableNames',false);
