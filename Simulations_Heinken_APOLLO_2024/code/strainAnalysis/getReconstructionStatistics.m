
% Get reconstruction statistics
mkdir([rootDir filesep 'data' filesep 'plot_ModelStatistics'])
cd([rootDir filesep 'data' filesep 'plot_ModelStatistics'])

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
    'Pasolli',[rootDir filesep 'data' filesep 'PasolliReconstructions']
    'Almeida',[rootDir filesep 'data' filesep 'AlmeidaReconstructions']
    };

versions = {
    '_draft','translatedDraftReconstructions'
    '_refined','refinedReconstructions'
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
            
            save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'All_statistics_' paths{p,1} versions{v,1} '.mat'],'stats');
            save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'UniqueReactions' paths{p,1} versions{v,1} '.mat'],'uniqueRxns')
            save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'UniqueMetabolites' paths{p,1} versions{v,1} '.mat'],'uniqueMets')
        end
        save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'All_statistics_' paths{p,1} versions{v,1} '.mat'],'stats');
        save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'UniqueReactions' paths{p,1} versions{v,1} '.mat'],'uniqueRxns')
        save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'UniqueMetabolites' paths{p,1} versions{v,1} '.mat'],'uniqueMets')
    end
end

%% add AGORA2 reconstructions

% To retrieve AGORA2 from the Virtual Metabolic Human website:
mkdir([rootDir filesep 'AGORA2'])
agora2_info = readInputTableForPipeline('AGORA2_infoFile.xlsx');
for i=2:size(agora2_info,1)
    websave([rootDir filesep 'AGORA2' filesep agora2_info{i,1} '.mat'],['https://www.vmh.life/files/reconstructions/AGORA2/version2.01/mat_files/individual_reconstructions/' agora2_info{i,1} '.mat'])
end

modelFolder=[rootDir filesep 'AGORA2'];
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

    save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'All_statistics_AGORA2.mat'],'stats');
    save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'UniqueReactions_AGORA2.mat'],'uniqueRxns')
    save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'UniqueMetabolites_AGORA2.mat'],'uniqueMets')
end
save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'All_statistics_AGORA2.mat'],'stats');
save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'UniqueReactions_AGORA2.mat'],'uniqueRxns')
save([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'UniqueMetabolites_AGORA2.mat'],'uniqueMets')
