% Get reconstruction statistics
rootDir=pwd;
mkdir([rootDir filesep 'Model_properties_analysis' filesep 'Overview_Model_stats'])
mkdir([rootDir filesep 'Model_properties_analysis' filesep 'Overview_Model_stats' filesep 'StatsByPhylum'])
cd([rootDir filesep 'Model_properties_analysis' filesep 'Overview_Model_stats'])

numWorkers=20;

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

global CBT_LP_SOLVER
solver = CBT_LP_SOLVER;

if numWorkers>0 && ~isempty(ver('parallel'))
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end
environment = getEnvironment();

% define diets
% Western diet
westernDiet = table2cell(readtable('WesternDietAGORA2.txt', 'Delimiter', 'tab'));
westernDiet=cellstr(string(westernDiet));

% High fiber diet
hfDiet = table2cell(readtable('HighFiberDietAGORA2.txt', 'Delimiter', 'tab'));
hfDiet=cellstr(string(hfDiet));

% define input parameters for findFluxConsistentSubset
param=struct;
feasTol = getCobraSolverParams('LP', 'feasTol');
param.('feasTol')= feasTol;
param.('epsilon')=feasTol*100;
param.('modeFlag')=0;
param.('method')='fastcc';

% define the paths to the two datasets
paths={
    '150k','C:\Users\Almut\OneDrive - National University of Ireland, Galway\150k_Project_Models'
    '90k','C:\Users\Almut\OneDrive - National University of Ireland, Galway\90000_genomes_Almeida_2019'
    };

versions = {
    '_draft','translatedDraftReconstructions'
    '_refined','curatedReconstructions'
    };

%% do computations
for p=1:size(paths,1)
    for v=1:size(versions,1)
        modelFolder = [paths{p,2} filesep versions{v,2}];
        dInfo = dir(modelFolder);
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
        stats{1,8}='Growth_aerobic_HFD';
        stats{1,9}='Growth_anaerobic_HFD';
        stats{1,10}='ATP_aerobic';
        stats{1,11}='ATP_anaerobic';
        stats{1,12}='Reactions';
        stats{1,13}='Metabolites';
        stats{1,14}='Genes';
        stats{1,15}='Flux consistent reactions';
        stats{1,16}='Stoich consistent reactions';

        % get unique reactions and metabolites
        uniqueRxns = {};
        uniqueMets = {};
        
        % proceed in batches for improved effiency
        steps=5000;
        for s=1:steps:length(modelList)
            if length(modelList)-s>=steps-1
                endPnt=steps-1;
            else
                endPnt=length(modelList)-s;
            end

            statsTmp={};
            uniqueRxnsTmp = {};
            uniqueMetsTmp = {};

            % Starting simulations
            parfor i=s:s+endPnt
                restoreEnvironment(environment);
                changeCobraSolver(solver, 'LP', 0, -1);

                modelIn=load([modelFolder filesep modelList{i}]);
                loadmodel=fieldnames(modelIn);
                model=modelIn.(loadmodel{1});

                % growth rates
                biomassID=find(strncmp(model.rxns,'bio',3));
                [AerobicGrowth, AnaerobicGrowth] = testGrowth(model, model.rxns(biomassID));
                statsTmp{i}{1}=strrep(modelList{i},'.mat','');
                statsTmp{i}{2}=AerobicGrowth(1,1);
                statsTmp{i}{3}=AnaerobicGrowth(1,1);
                statsTmp{i}{4}=AerobicGrowth(1,2);
                statsTmp{i}{5}=AnaerobicGrowth(1,2);

                % Western diet
                model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), -1000, 'l');
                model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), 1000, 'u');
                model = useDiet(model,westernDiet);
                model = changeRxnBounds(model, 'EX_o2(e)', -10, 'l');
                FBA=optimizeCbModel(model,'max');
                statsTmp{i}{6}=FBA.f;
                model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
                FBA=optimizeCbModel(model,'max');
                statsTmp{i}{7}=FBA.f;

                % High fiber diet
                model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), -1000, 'l');
                model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), 1000, 'u');
                model = useDiet(model,hfDiet);
                model = changeRxnBounds(model, 'EX_o2(e)', -10, 'l');
                FBA=optimizeCbModel(model,'max');
                statsTmp{i}{8}=FBA.f;
                model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
                FBA=optimizeCbModel(model,'max');
                statsTmp{i}{9}=FBA.f;

                % ATP
                [ATPFluxAerobic, ATPFluxAnaerobic] = testATP(model);
                statsTmp{i}{10}=ATPFluxAerobic(1,1);
                statsTmp{i}{11}=ATPFluxAnaerobic(1,1);

                % Number of reactions, metabolites, and genes
                statsTmp{i}{12}=length(model.rxns);
                statsTmp{i}{13}=length(model.mets);
                statsTmp{i}{14}=length(model.genes);

                % stoichiometric and flux consistency
                [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool] = findFluxConsistentSubset(model,param);
                % exclude exchange and demand reactions
                exRxns=vertcat(find(strncmp(model.rxns,'EX_',3)),find(strcmp(model.rxns,'rxn00062')));
                fluxConsistentRxnBool(exRxns,:)=[];
                fluxInConsistentRxnBool(exRxns,:)=[];
                statsTmp{i}{15}=sum(fluxConsistentRxnBool)/(sum(fluxConsistentRxnBool) + sum(fluxInConsistentRxnBool));
                [SConsistentMetBool,SConsistentRxnBool,SInConsistentMetBool,SInConsistentRxnBool,unknownSConsistencyMetBool,unknownSConsistencyRxnBool]=...
                    findStoichConsistentSubset(model);
                % exclude exchange and demand reactions
                SConsistentRxnBool(exRxns,:)=[];
                SInConsistentRxnBool(exRxns,:)=[];
                statsTmp{i}{16}=sum(SConsistentRxnBool)/(sum(SConsistentRxnBool) + sum(SInConsistentRxnBool));
                
                uniqueRxnsTmp = union(uniqueRxnsTmp,model.rxns);
                uniqueMetsTmp = union(uniqueMetsTmp,model.mets);
            end

            for i=s:s+endPnt
                for n=1:16
                    stats{i+1,n}=statsTmp{i}{n};
                end
            end

            uniqueRxns = union(uniqueRxns,uniqueRxnsTmp);
            uniqueMets = union(uniqueMets,uniqueMetsTmp);

            save(['All_statistics_' paths{p,1} versions{v,1} '.mat'],'stats');
            save(['UniqueReactions' paths{p,1} versions{v,1} '.mat'],'uniqueRxns')
            save(['UniqueMetabolites' paths{p,1} versions{v,1} '.mat'],'uniqueMets')
        end
        cell2csv(['All_statistics_' paths{p,1} versions{v,1} '.csv'],stats);
        save(['UniqueReactions' paths{p,1} versions{v,1} '.mat'],'uniqueRxns')
        save(['UniqueMetabolites' paths{p,1} versions{v,1} '.mat'],'uniqueMets')
    end
end

%% export some statistics on taxa in the resource
% get statistics for taxon presence
taxa={'Phylum','Class','Order','Family','Genus'};
taxonStats={'','Pasolli','Almeida','GlobalRecon'};
taxonStats{2,1}='Strains';
for i=1:length(taxa)
    taxonStats{i+2,1}=taxa{i};
end
taxonStats{8,1}='Named species';

for p=1:size(paths,1)
taxonomy = readInputTableForPipeline([paths{p,2} filesep 'inputFiles' filesep 'adaptedInfoFile.txt']);
    taxonStats{2,p+1}=length(unique(taxonomy(2:end,1)));
    for i=1:length(taxa)
        taxCol=find(strcmp(taxonomy(1,:),taxa{i}));
        getTax=unique(taxonomy(2:end,taxCol));
        getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
        taxonStats{i+2,p+1}=length(getTax);
    end
    % get number of named species
    taxCol=find(strcmp(taxonomy(1,:),'Species'));
    getTax=unique(taxonomy(2:end,taxCol));
    getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
    chrTax=getTax;
    chrTax(find(strcmp(chrTax,'')),:)=[];
    chrTax(find(contains(chrTax,'sp.')),:)=[];
    chrTax(find(contains(chrTax,'nov.')),:)=[];
    chrTax(find(contains(chrTax,'cf')),:)=[];
    chrTax(find(contains(chrTax,'iales')),:)=[];
    chrTax(find(contains(chrTax,'aceae')),:)=[];
    chrTax(find(contains(chrTax,'uncultured')),:)=[];
    chrTax(find(contains(chrTax,'taxon')),:)=[];
    chrTax(find(strncmp(chrTax,'bacterium',9)),:)=[];
    chrTax(find(contains(chrTax,' sp')),:)=[];
    taxonStats{8,p+1}=length(chrTax);

    %% separate the results by phylum
    % define path to the taxonomy information
    taxonomy = readInputTableForPipeline([paths{p,2} filesep 'inputFiles' filesep 'adaptedInfoFile.txt']);
    [phyla,~,J]=unique(taxonomy(2:end,find(strcmp(taxonomy(1,:),'Phylum'))));
    % remove the phyla with only few entries and unclassified bacteria
    cnt = histc(J, 1:numel(phyla));
    phyla(cnt<3)=[];
    phyla(strncmp(phyla,'unclassified',length('unclassified')))=[];
    phyla(strncmp(phyla,'Unclassified',length('Unclassified')))=[];
    phyla(strncmp(phyla,'N/A',3))=[];
    phyla(strcmp(phyla,''))=[];

    % get the statistics by phylum
    stats = readInputTableForPipeline(['All_statistics_' paths{p,1} '_refined.csv']);
    for i=2:size(stats,2)
        statsByPhylum={};
        for j=1:length(phyla)
            statsByPhylum{j,1}=phyla{j};
            allbacs=taxonomy(find(strcmp(taxonomy(:,find(strcmp(taxonomy(1,:),'Phylum'))),phyla{j})),1);
            for k=1:length(allbacs)
                statsByPhylum{j,k+1}=stats{find(strcmp(stats(:,1),allbacs{k})),i};
            end
        end
        statsByPhylum=flip(statsByPhylum);
        statsByPhylum=statsByPhylum';
        cell2csv(['StatsByPhylum' filesep stats{1,i} '_' paths{p,1} '.csv'],statsByPhylum);
    end
end

% get taxon stats for all reconstructions taken together
taxonomyAll = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Combined_taxonomy_info.xlsx']);
taxonStats{2,4}=length(unique(taxonomyAll(2:end,1)));
for i=1:length(taxa)
    taxCol=find(strcmp(taxonomyAll(1,:),taxa{i}));
    getTax=unique(taxonomyAll(2:end,taxCol));
    getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
    taxonStats{i+2,4}=length(getTax);
end
% get number of characterized and uncharacterized species
taxCol=find(strcmp(taxonomyAll(1,:),'Species'));
getTax=unique(taxonomyAll(2:end,taxCol));
getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
chrTax=getTax;
chrTax(find(contains(chrTax,'sp.')),:)=[];
chrTax(find(contains(chrTax,'nov.')),:)=[];
chrTax(find(contains(chrTax,'cf')),:)=[];
chrTax(find(contains(chrTax,'iales')),:)=[];
chrTax(find(contains(chrTax,'aceae')),:)=[];
chrTax(find(contains(chrTax,'uncultured')),:)=[];
chrTax(find(contains(chrTax,'taxon')),:)=[];
chrTax(find(strncmp(chrTax,'bacterium',9)),:)=[];

taxonStats{8,4}=length(chrTax);
taxonStats{9,4}=length(getTax) - length(chrTax);
cell2csv('Taxon_statistics.csv',taxonStats);

%% add AGORA2 reconstructions
modelFolder='C:\Users\Almut\National University of Ireland, Galway\Group_MSP - Documents\AGORA2\Current_Version_AGORA2\Output_Models';
dInfo = dir(modelFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(find(~contains(modelList,'.mat')),:)=[];

stats={};
stats{1,1}='Model_ID';
stats{1,2}='Growth_aerobic_UM';
stats{1,3}='Growth_anaerobic_UM';
stats{1,4}='Growth_aerobic_CM';
stats{1,5}='Growth_anaerobic_CM';
stats{1,6}='Growth_aerobic_WD';
stats{1,7}='Growth_anaerobic_WD';
stats{1,8}='Growth_aerobic_HFD';
stats{1,9}='Growth_anaerobic_HFD';
stats{1,10}='ATP_aerobic';
stats{1,11}='ATP_anaerobic';
stats{1,12}='Reactions';
stats{1,13}='Metabolites';
stats{1,14}='Genes';
stats{1,15}='Flux consistent reactions';
stats{1,16}='Stoich consistent reactions';

% get unique reactions and metabolites
uniqueRxns = {};
uniqueMets = {};

% proceed in batches for improved effiency
steps=1000;
for s=1:steps:length(modelList)
    if length(modelList)-s>=steps-1
        endPnt=steps-1;
    else
        endPnt=length(modelList)-s;
    end

    statsTmp={};
    uniqueRxnsTmp = {};
    uniqueMetsTmp = {};

    % Starting simulations
    parfor i=s:s+endPnt
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);

        modelIn=load([modelFolder filesep modelList{i}]);
        loadmodel=fieldnames(modelIn);
        model=modelIn.(loadmodel{1});

        % growth rates
        biomassID=find(strncmp(model.rxns,'bio',3));
        [AerobicGrowth, AnaerobicGrowth] = testGrowth(model, model.rxns(biomassID));
        statsTmp{i}{1}=strrep(modelList{i},'.mat','');
        statsTmp{i}{2}=AerobicGrowth(1,1);
        statsTmp{i}{3}=AnaerobicGrowth(1,1);
        statsTmp{i}{4}=AerobicGrowth(1,2);
        statsTmp{i}{5}=AnaerobicGrowth(1,2);

        % Western diet
        model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), -1000, 'l');
        model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), 1000, 'u');
        model = useDiet(model,westernDiet);
        model = changeRxnBounds(model, 'EX_o2(e)', -10, 'l');
        FBA=optimizeCbModel(model,'max');
        statsTmp{i}{6}=FBA.f;
        model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
        FBA=optimizeCbModel(model,'max');
        statsTmp{i}{7}=FBA.f;

        % High fiber diet
        model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), -1000, 'l');
        model = changeRxnBounds(model, model.rxns(strncmp('EX_', model.rxns, 3)), 1000, 'u');
        model = useDiet(model,hfDiet);
        model = changeRxnBounds(model, 'EX_o2(e)', -10, 'l');
        FBA=optimizeCbModel(model,'max');
        statsTmp{i}{8}=FBA.f;
        model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
        FBA=optimizeCbModel(model,'max');
        statsTmp{i}{9}=FBA.f;

        % ATP
        [ATPFluxAerobic, ATPFluxAnaerobic] = testATP(model);
        statsTmp{i}{10}=ATPFluxAerobic(1,1);
        statsTmp{i}{11}=ATPFluxAnaerobic(1,1);

        % Number of reactions, metabolites, and genes
        statsTmp{i}{12}=length(model.rxns);
        statsTmp{i}{13}=length(model.mets);
        statsTmp{i}{14}=length(model.genes);

        % stoichiometric and flux consistency
        [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool] = findFluxConsistentSubset(model,param);
        % exclude exchange and demand reactions
        exRxns=vertcat(find(strncmp(model.rxns,'EX_',3)),find(strcmp(model.rxns,'rxn00062')));
        fluxConsistentRxnBool(exRxns,:)=[];
        fluxInConsistentRxnBool(exRxns,:)=[];
        statsTmp{i}{15}=sum(fluxConsistentRxnBool)/(sum(fluxConsistentRxnBool) + sum(fluxInConsistentRxnBool));
        [SConsistentMetBool,SConsistentRxnBool,SInConsistentMetBool,SInConsistentRxnBool,unknownSConsistencyMetBool,unknownSConsistencyRxnBool]=...
            findStoichConsistentSubset(model);
        % exclude exchange and demand reactions
        SConsistentRxnBool(exRxns,:)=[];
        SInConsistentRxnBool(exRxns,:)=[];
        statsTmp{i}{16}=sum(SConsistentRxnBool)/(sum(SConsistentRxnBool) + sum(SInConsistentRxnBool));

        uniqueRxnsTmp = union(uniqueRxnsTmp,model.rxns);
        uniqueMetsTmp = union(uniqueMetsTmp,model.mets);
    end

    for i=s:s+endPnt
        for n=1:16
            stats{i+1,n}=statsTmp{i}{n};
        end
    end
    uniqueRxns = union(uniqueRxns,uniqueRxnsTmp);
    uniqueMets = union(uniqueMets,uniqueMetsTmp);

    save('All_statistics_AGORA2.mat','stats');
    save('UniqueReactions_AGORA2.mat','uniqueRxns')
    save('UniqueMetabolites_AGORA2.mat','uniqueMets')
end
cell2csv('All_statistics_AGORA2.csv',stats);
save('UniqueReactions_AGORA2.mat','uniqueRxns')
save('UniqueMetabolites_AGORA2.mat','uniqueMets')

%% create a table with stats on both resources separately and together
table={'','Almeida','Pasolli','GlobalBiome','AGORA2'
    'Number of reconstructions','','','',''
    'Aerobic growth','','','',''
    'Anaerobic growth','','','',''
    'Aerobic ATP production','','','',''
    'Anaerobic ATP production','','','',''
    'Number of reactions','','','',''
    'Number of metabolites','','','',''
    'Number of genes','','','',''
    'Flux consistent reactions','','','',''
    'Stoich consistent reactions','','','',''
    'Unique reactions','','','',''
    'Unique metabolites','','','',''
    };

% load data
data_150k = readInputTableForPipeline([rootDir filesep 'Overview_Model_stats' filesep 'All_statistics_' paths{1,1} '_refined.csv']);
data_150k(1,:) = [];
data_90k = readInputTableForPipeline([rootDir filesep 'Overview_Model_stats' filesep 'All_statistics_' paths{2,1} '_refined.csv']);
data_90k(1,:) = [];
data_combined=vertcat(data_150k,data_90k(1:end,:));
data_AGORA2 = readInputTableForPipeline([rootDir filesep 'Overview_Model_stats' filesep 'All_statistics_AGORA2.csv']);
data_AGORA2(1,:) = [];
data_150k = cell2mat(data_150k(:,2:end));
data_90k = cell2mat(data_90k(:,2:end));
data_combined = cell2mat(data_combined(:,2:end));
data_AGORA2 = cell2mat(data_AGORA2(:,2:end));

table{2,2}=length(data_90k);
table{2,3}=length(data_150k);
table{2,4}=length(data_combined);
table{2,5}=length(data_AGORA2);

% there may be NaNs in the data
data_90k(find(isnan(data_90k(:,6))),:)=[];
data_90k(find(isnan(data_90k(:,7))),:)=[];
data_150k(find(isnan(data_150k(:,6))),:)=[];
data_150k(find(isnan(data_150k(:,7))),:)=[];
data_combined(find(isnan(data_combined(:,6))),:)=[];
data_combined(find(isnan(data_combined(:,7))),:)=[];

for i=3:4
% 90k averages
av = mean(data_90k(:,i+1));
s = std(data_90k(:,i+1));
table{i,2} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
% 150k averages
av = mean(data_150k(:,i+1));
s = std(data_150k(:,i+1));
table{i,3} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
% total averages
av = mean(data_combined(:,i+1));
s = std(data_combined(:,i+1));
table{i,4} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
% AGORA2 averages
av = mean(data_AGORA2(:,i+1));
s = std(data_AGORA2(:,i+1));
table{i,5} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
end

for i=5:11
    % 90k averages
    av = mean(data_90k(:,i+4));
    s = std(data_90k(:,i+4));
    table{i,2} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
    % 150k averages
    av = mean(data_150k(:,i+4));
    s = std(data_150k(:,i+4));
    table{i,3} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
    % total averages
    av = mean(data_combined(:,i+4));
    s = std(data_combined(:,i+4));
    table{i,4} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
    % AGORA2 averages
    av = mean(data_AGORA2(:,i+4));
    s = std(data_AGORA2(:,i+4));
    table{i,5} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
end
cell2csv('Summary_ReconStats.csv',table);

%% get stats on taxon distribution

taxTable={'','Almeida','Pasolli','GlobalBiome','AGORA2'
    'Phyla','','','',''
    'Classes','','','',''
    'Order','','','',''
    'Families','','','',''
    'Genera','','','',''
    'Species','','','',''
    'Unclassified on species level','','','',''
};
taxa={'Phylum','Class','Order','Family','Genus','Species'};

infoFiles = {
    [rootDir filesep 'Model_properties_analysis' filesep 'inputFiles' filesep 'mags-gut_qs50_checkm_adaptedToPipeline.txt']
    [rootDir filesep 'Model_properties_analysis' filesep 'inputFiles' filesep 'SequencedGenomesInfo.txt']
    [rootDir filesep 'Model_properties_analysis' filesep 'inputFiles' filesep 'Combined_taxonomy_info.txt']
    'AGORA2_infoFile.xlsx'
    };

for i=1:length(infoFiles)
    info = readInputTableForPipeline(infoFiles{i});
    for j=1:length(taxa)
        taxCol=find(strcmp(info(1,:),taxa{j}));
        getTax=unique(info(2:end,taxCol));
        getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
        getTax(find(strcmp(getTax,'')),:)=[];
        getTax(find(strcmp(getTax,'N/A')),:)=[];
        getTax(find(strcmp(getTax,'NA')),:)=[];
        taxTable{j+1,i+1}=length(getTax);
    end
    % get unclassified at least on species level
    taxCol=find(strcmp(info(1,:),taxa{j}));
    uncl=0;
    getUncl=find(strncmp(info(2:end,taxCol),'unclassified',length('unclassified')));
    uncl = uncl + length(getUncl);
    getUncl=find(strcmp(info(2:end,taxCol),''));
    uncl = uncl + length(getUncl);
    getUncl=find(strcmp(info(2:end,taxCol),'N/A'));
    uncl = uncl + length(getUncl);
    getUncl=find(strcmp(info(2:end,taxCol),'NA'));
    uncl = uncl + length(getUncl);
    taxTable{8,i+1} = uncl/(size(info,1)-1);
end
cell2csv('Summary_TaxStats.csv',table);

%% separate the results by phylum
taxonomyAll = readInputTableForPipeline([rootDir filesep 'inputFiles' filesep 'Combined_taxonomy_info.xlsx']);

[phyla,~,J]=unique(taxonomyAll(2:end,find(strcmp(taxonomyAll(1,:),'Phylum'))));
% remove the phyla with only few entries and unclassified bacteria
cnt = histc(J, 1:numel(phyla));
phyla(cnt<3)=[];
phyla(strncmp(phyla,'unclassified',length('unclassified')))=[];
phyla(strncmp(phyla,'Unclassified',length('Unclassified')))=[];
phyla(strncmp(phyla,'N/A',3))=[];
phyla(strcmp(phyla,''))=[];

% get the statistics by phylum
% load data
data_150k = readInputTableForPipeline(['All_statistics_' paths{1,1} '_refined.csv']);
data_90k = readInputTableForPipeline(['All_statistics_' paths{2,1} '_refined.csv']);
data_combined=vertcat(data_150k,data_90k(2:end,:));

for i=2:size(data_combined,2)
    statsByPhylum={};
    for j=1:length(phyla)
        statsByPhylum{j,1}=phyla{j};
        allbacs=taxonomyAll(find(strcmp(taxonomyAll(:,find(strcmp(taxonomyAll(1,:),'Phylum'))),phyla{j})),1);
        for k=1:length(allbacs)
            statsByPhylum{j,k+1}=data_combined{find(strcmp(data_combined(:,1),allbacs{k})),i};
        end
    end
    statsByPhylum=flip(statsByPhylum);
    statsByPhylum=statsByPhylum';
    cell2csv(['StatsByPhylum' filesep data_combined{1,i} '_combined.csv'],statsByPhylum);
end
