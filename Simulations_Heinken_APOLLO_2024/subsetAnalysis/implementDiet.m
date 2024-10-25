
% create microbiome models with diet constraints

numWorkers = 12;
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

%%
% modPath = [rootDir filesep 'MicrobiomeModeling' filesep 'MicrobiomeModels'];
modPath = [pwd filesep 'MicrobiomeModels'];

dInfo = dir(modPath);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'microbiota_model'),:)=[];

dietFilePath='AverageEuropeanDiet_modified';
HumanMets={'gchola','-10';'tdchola','-10';'tchola','-10';'dgchol','-10';'34dhphe','-10';'5htrp','-10';'Lkynr','-10';'f1a','-1';'gncore1','-1';'gncore2','-1';'dsT_antigen','-1';'sTn_antigen','-1';'core8','-1';'core7','-1';'core5','-1';'core4','-1';'ha','-1';'cspg_a','-1';'cspg_b','-1';'cspg_c','-1';'cspg_d','-1';'cspg_e','-1';'hspg','-1'};
lowerBMBound=0;
upperBMBound=1;
% dietPath = [rootDir filesep 'MicrobiomeModeling' filesep 'Microbiome_models_AED'];
dietPath = [pwd filesep 'Microbiome_models_AED'];

mkdir(dietPath)

uGrowth = {};
dietGrowth = {};

%%

[diet] = adaptVMHDietToAGORA(dietFilePath,'Microbiota');

% proceed in batches for improved effiency
steps=50;
for s=1:steps:length(modelList)
    if length(modelList)-s>=steps-1
        endPnt=steps-1;
    else
        endPnt=length(modelList)-s;
    end
    
    presolTmp={};
    infesMatTmp={};
    
    % Starting personalized simulations
    parfor k=s:s+endPnt
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);
        
        sampleID = modelList{k,1};
        
        % microbiota_model=readCbModel[modPath filesep (modelList{k}]);
        modelStr=load([modPath filesep modelList{k}]);
        modelF=fieldnames(modelStr);
        model=modelStr.(modelF{1});

        for j = 1:length(model.rxns)
            if strfind(model.rxns{j}, 'biomass')
                model.lb(j) = 0;
            end
        end

        % adapt constraints
        BiomassNumber=find(strcmp(model.rxns,'communityBiomass'));
        Components = model.mets(find(model.S(:, BiomassNumber)));
        Components = strrep(Components,'_PGP[c]','');
        for j=1:length(Components)
            % remove constraints on demand reactions to prevent infeasibilities
            findDm= model.rxns(find(strncmp(model.rxns,[Components{j} '_DM_'],length([Components{j} '_DM_']))));
            model = changeRxnBounds(model, findDm, 0, 'l');
            % constrain flux through sink reactions
            findSink= model.rxns(find(strncmp(model.rxns,[Components{j} '_sink_'],length([Components{j} '_sink_']))));
            model = changeRxnBounds(model, findSink, -1, 'l');
        end

        model = changeObjective(model, 'EX_microbeBiomass[fe]');
        AllRxn = model.rxns;
        RxnInd = find(cellfun(@(x) ~isempty(strfind(x, '[d]')), AllRxn));
        EXrxn = model.rxns(RxnInd);
        EXrxn = regexprep(EXrxn, 'EX_', 'Diet_EX_');
        model.rxns(RxnInd) = EXrxn;
        model = changeRxnBounds(model, 'communityBiomass', lowerBMBound, 'l');
        model = changeRxnBounds(model, 'communityBiomass', upperBMBound, 'u');
        model=changeRxnBounds(model,model.rxns(strmatch('UFEt_',model.rxns)),1000000,'u');
        model=changeRxnBounds(model,model.rxns(strmatch('DUt_',model.rxns)),1000000,'u');
        model=changeRxnBounds(model,model.rxns(strmatch('EX_',model.rxns)),1000000,'u');

        % Using input diet
        model_sd=model;

        [model_sd] = useDiet(model_sd, diet,0);

        % add the human metabolites
        for l=1:length(HumanMets)
            model_sd=changeRxnBounds(model_sd,strcat('Diet_EX_',HumanMets{l},'[d]'),str2num(HumanMets{l,2}),'l');
        end

        % test growth
        solution_uDiet=optimizeCbModel(model);
        solution_sDiet=optimizeCbModel(model_sd);
        uGrowthTmp{k}=solution_uDiet.f;
        dietGrowthTmp{k}=solution_sDiet.f;

        % apply lower bound if possible
        if solution_sDiet.f >0.4
            model_sd=changeRxnBounds(model_sd,'communityBiomass',0.4,'l');
        end

        microbiota_model=model_sd;
        parsave([dietPath filesep modelList{k}],microbiota_model)
    end
    
    for k=s:s+endPnt
        uGrowth{k,1} = modelList{k};
        uGrowth{k,2} = uGrowthTmp{k};
        dietGrowth{k,1} = modelList{k};
        dietGrowth{k,2} = dietGrowthTmp{k};
    end
    save('uGrowth.mat','uGrowth')
    save('dietGrowth.mat','dietGrowth')
end
