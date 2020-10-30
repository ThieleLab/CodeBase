%%
%script to extract flux bounds from gene expression data based on
%differential gene expression study
%%
%study A // 451 genes
%http://www.ncbi.nlm.nih.gov/pubmed/19912253
geneList = readtable('diffGenesPancreasBiodbJp.csv');
entrezList = table2array(geneList(:,2));
%%
%load files
harvey = modelOrganAllCoupled;
changeCobraSolver('tomlab_cplex');
%%
entrezList = entrezList(~isnan(entrezList));
entrezList = num2cell(entrezList);
for i = 1:length(entrezList)
    entrezList{i} = num2str(entrezList{i});
end
for i = 1:length(harvey.genes)
     nGene = strsplit(harvey.genes{i},'.');
     harveyGenes{i} = char(nGene(1));
end
[studyaGenes,studyaInd] = intersect(entrezList,harveyGenes');
geneList(studyaInd,1)
%regvec is the gene fold change value for the selected genes (extracted from
%the paper in regLevels.txt)
regVec = [-1.58 -1.64 -0.6075 -3.19 1.03 -1.55 -4.72 -2.195  1.42...
    1.4 2.28 1.525 -2.565 -1.08 1.95 1.43 1.686 1.44 1.66 1.383 -4.5...
    1.54 -2.28 2.025];
fprintf('There are %d metabolic genes among %d differentially expressed genes ',length(regVec),length(entrezList));
%%
%load healthy patients fluxes (take FVA solution)
load fvaminHealthyrev;%fvaminH
load fvamaxHealthyrev;%fvamaxH
%deriving constraints
%genes involved in T1D
consGenes = [433;5;8;147;371;11;178;14;16;266;355;349;148;338;343;3;...
    348;446;177;455;1;454;342;78];
%fold change of genes
geneExpFoldChange = [-1.58 -1.64 -0.6075 -3.19 1.03 -1.55 -4.72 -2.195  1.42...
1.4 2.28 1.525 -2.565 -1.08 1.95 1.43 1.686 1.44 1.66 1.383 -4.5...
1.54 -2.28 2.025];
geneExpFoldChangeQuan = 2.^geneExpFoldChange;
rxnFlux = [];
geneCorrespondance = [];
for i=consGenes'
    rxnFlux = [rxnFlux;find(harvey.rxnGeneMat(:,i))];
    geneCorrespondance = [geneCorrespondance;repmat(i,length(find(harvey.rxnGeneMat(:,i))),1)];
end
allRxns = length(rxnFlux);
[rxnFlux,fluxIndex] = unique(rxnFlux);
geneCorrespondance = geneCorrespondance(fluxIndex);
fprintf('There are %d reactions that corrsepond to the genes and only %d are unique \n', allRxns ,length(rxnFlux));
%here we average the fluxes in the healthy state and derive constraints
%get average reaction flux for the healthy individuals
sum = mean([fvaminH(rxnFlux,:) fvamaxH(rxnFlux,:)],2);%2076
%see intersection with active reactions
fprintf('There are %d active reactions among %d possible reactions ', length(find(sum)) ,length(sum));
%compute new constraints
newConstraints = zeros(length(find(sum)),1);
fc = zeros(length(find(sum)),1);
for i =1:length(find(sum))
    c = find(consGenes == geneCorrespondance(i));
    newConstraints(i) = sum(i)*geneExpFoldChangeQuan(c);
    fc(i)=geneExpFoldChangeQuan(c);
end
%take only pancreas fluxes
pancreasId = cellfun(@(x) ~isempty(strmatch(x,'Pancreas')),harvey.organs(rxnFlux));
GeneExCons.rxnIdPancreas = rxnFlux(pancreasId);
GeneExCons.newConstraintsPancreas = newConstraints(pancreasId);
GeneExCons.deletedCons=[];
GeneExCons.fc=fc(pancreasId);
fprintf('\n There are %d reactions that take place in the pancreas ', length(GeneExCons.rxnIdPancreas));
cd(['..' filesep '..' filesep 'data'])
%save data file
save('GeneExCons.mat','GeneExCons');