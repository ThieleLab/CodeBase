% Get reconstruction statistics for Figure 1

mkdir('Data_Figure_1')

% Get the taxonomy information
taxonomy = readInputTableForPipeline('AGORA2_infoFile.xlsx');

%% Curated reconstructions
resultsFolder=[rootDir filesep 'Data_Figure_1'];
dInfo = dir(refinedFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

% Western diet
westernDiet = table2cell(readtable('WesternDietAGORA2.txt', 'Delimiter', 'tab'));
westernDiet=cellstr(string(westernDiet));

% get reconstruction statistics: stats, stats, production, number of
% reactions, metabolites, and genes

stats={};
stats{1,1}='Model_ID';
stats{1,2}='Growth_aerobic_UM';
stats{1,3}='Growth_anaerobic_UM';
stats{1,4}='Growth_aerobic_CM';
stats{1,5}='Growth_anaerobic_CM';
stats{1,6}='Growth_aerobic_WD';
stats{1,7}='Growth_anaerobic_WD';
stats{1,8}='ATP_aerobic';
stats{1,9}='ATP_anaerobic';
stats{1,10}='Reactions';
stats{1,11}='Metabolites';
stats{1,12}='Genes';

for i=1:length(modelList)
    i
    model=readCbModel([refinedFolder filesep modelList{i}]);
    % growth rates
    biomassID=find(strncmp(model.rxns,'bio',3));
    [AerobicGrowth, AnaerobicGrowth] = testGrowth(model, model.rxns(biomassID));
    stats{i+1,1}=strrep(modelList{i},'.mat','');
    stats{i+1,2}=AerobicGrowth(1,1);
    stats{i+1,3}=AnaerobicGrowth(1,1);
    stats{i+1,4}=AerobicGrowth(1,2);
    stats{i+1,5}=AnaerobicGrowth(1,2);
    % Western diet
    model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), -1000, 'l');
    model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), 1000, 'u');
    model = useDiet(model,westernDiet);
    model = changeRxnBounds(model, 'EX_o2(e)', -10, 'l');
    FBA=optimizeCbModel(model,'max');
    stats{i+1,6}=FBA.f;
    model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
    FBA=optimizeCbModel(model,'max');
    stats{i+1,7}=FBA.f;
    % ATP
    [ATPFluxAerobic, ATPFluxAnaerobic] = testATP(model);
    stats{i+1,8}=ATPFluxAerobic(1,1);
    stats{i+1,9}=ATPFluxAnaerobic(1,1);
    % Number of reactions, metabolites, and genes
    stats{i+1,10}=length(model.rxns);
    stats{i+1,11}=length(model.mets);
    stats{i+1,12}=length(model.genes);
end
cell2csv([resultsFolder filesep 'All_statistics_Curated.csv'],stats);

%% separate the results by phylum

[phyla,~,J]=unique(taxonomy(2:end,find(strcmp(taxonomy(1,:),'Phylum'))));
% remove the phyla with only few entries and unclassified bacteria
cnt = histc(J, 1:numel(phyla));
phyla(cnt<3)=[];
phyla(strncmp(phyla,'unclassified',length('unclassified')))=[];
phyla(strncmp(phyla,'Unclassified',length('Unclassified')))=[];

% get the statistics by phylum
for i=2:size(stats,2)
    statsByPhylum={};
    for j=1:length(phyla)
        statsByPhylum{j,1}=phyla{j};
        allbacs=taxonomy(find(strcmp(taxonomy(:,find(strcmp(taxonomy(1,:),'Phylum'))),phyla{j})),1);
        for k=1:length(allbacs)
            statsByPhylum{j,k+1}=stats{find(strcmp(stats(:,1),allbacs{k})),i};
        end
    end
    % get remaining phyla and summarize them
    otherphyla=setdiff(taxonomy(2:end,find(strcmp(taxonomy(1,:),'Phylum'))),phyla);
    otherbacs=taxonomy(find(ismember(taxonomy(:,find(strcmp(taxonomy(1,:),'Phylum'))),otherphyla)),1);
    cnt=size(statsByPhylum,1)+1;
    statsByPhylum{cnt,1}='Other';
    for k=1:length(otherbacs)
        statsByPhylum{cnt,k+1}=stats{find(strcmp(stats(:,1),otherbacs{k})),i};
    end
    statsByPhylum=flip(statsByPhylum);
    statsByPhylum=statsByPhylum';
    cell2csv([resultsFolder filesep stats{1,i} '.csv'],statsByPhylum);
end

%% Draft reconstructions
dInfo = dir(translatedDraftsFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

% get reconstruction statistics: stats, stats, production, number of
% reactions, metabolites, and genes

stats={};
stats{1,1}='Model_ID';
stats{1,2}='Growth_aerobic_UM';
stats{1,3}='Growth_anaerobic_UM';
stats{1,4}='Growth_aerobic_CM';
stats{1,5}='Growth_anaerobic_CM';
stats{1,6}='Growth_aerobic_WD';
stats{1,7}='Growth_anaerobic_WD';
stats{1,8}='ATP_aerobic';
stats{1,9}='ATP_anaerobic';
stats{1,10}='Reactions';
stats{1,11}='Metabolites';
stats{1,12}='Genes';

for i=1:length(modelList)
    i
    model=readCbModel([translatedDraftsFolder filesep modelList{i}]);
    % growth rates
    biomassID=find(strncmp(model.rxns,'bio',3));
    [AerobicGrowth, AnaerobicGrowth] = testGrowth(model, model.rxns(biomassID));
    stats{i+1,1}=strrep(modelList{i},'.mat','');
    stats{i+1,2}=AerobicGrowth(1,1);
    stats{i+1,3}=AnaerobicGrowth(1,1);
    stats{i+1,4}=AerobicGrowth(1,2);
    stats{i+1,5}=AnaerobicGrowth(1,2);
    % Western diet
    model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), -1000, 'l');
    model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), 1000, 'u');
    model = useDiet(model,westernDiet);
    model = changeRxnBounds(model, 'EX_o2(e)', -10, 'l');
    FBA=optimizeCbModel(model,'max');
    stats{i+1,6}=FBA.f;
    model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
    FBA=optimizeCbModel(model,'max');
    stats{i+1,7}=FBA.f;
    % ATP
    [ATPFluxAerobic, ATPFluxAnaerobic] = testATP(model);
    stats{i+1,8}=ATPFluxAerobic(1,1);
    stats{i+1,9}=ATPFluxAnaerobic(1,1);
    % Number of reactions, metabolites, and genes
    stats{i+1,10}=length(model.rxns);
    stats{i+1,11}=length(model.mets);
    stats{i+1,12}=length(model.genes);
end
cell2csv([resultsFolder filesep 'All_statistics_Draft.csv'],stats);
