
%% Compute drug yields

resultsFolder = [rootDir filesep 'ComputedDrugFluxes' filesep];

database=loadVMHDatabase;

% get organisms with comparative genomics
genomeAnnotation = readtable([inputDataFolder filesep 'gapfilledGenomeAnnotation.txt'], 'ReadVariableNames', false, 'Delimiter', 'tab','TreatAsEmpty',['UND. -60001','UND. -2011','UND. -62011']);
genomeAnnotation = table2cell(genomeAnnotation);
orgs=unique(genomeAnnotation(:,1));

global CBT_LP_SOLVER
solver = CBT_LP_SOLVER;
% initialize parallel pool
if numWorkers > 0
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end
environment = getEnvironment();

%% compute yields

yieldsToCompute={
    'DM_atp_c_'
    'DM_co2[c]'
    'DM_pyr[c]'
    'DM_nh4[c]'
    };

basicCompounds={
    'EX_h2o(e)' '-100'
    'EX_h(e)' '-100'
    'EX_pi(e)' '-10'
    'EX_o2(e)' '-10'
    };

drugExchanges={'glc_D','EX_glc_D(e)';'sn38g','EX_sn38g(e)';'dfdcytd','EX_dfdcytd(e)';'fcsn','EX_fcsn(e)';'5fura','EX_5fura(e)';'r788','EX_r788(e)';'bzd','EX_bzd(e)';'lactl','EX_lactl(e)';'chlphncl','EX_chlphncl(e)';'5asa','EX_5asa(e)';'digoxin','EX_digoxin(e)';'srv','EX_srv(e)';'34dhphe','EX_34dhphe(e)';'tchola','EX_tchola(e)';'dopa','EX_dopa(e)'};

for k=1:length(yieldsToCompute)
    met=strrep(yieldsToCompute{k},'DM_','');
    met=strrep(met,'sink_','');
    met=strrep(met,'[c]','');
    met=strrep(met,'_c_','');
    
    fluxesTmp={};
    parfor i = 1:length(orgs)
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP');
        % prevent creation of log files
        changeCobraSolverParams('LP', 'logFile', 0);
        
        % load and constrain the model
        model=readCbModel([refinedFolder filesep orgs{i} '.mat']);
        model=addDemandReaction(model,{'co2[c]','pyr[c]','nh4[c]'});
        model = useDiet(model,basicCompounds);
        model.lb(find(strncmp(model.rxns,'sink_',5)))=0;
        model=changeObjective(model,yieldsToCompute{k});
        
        % compute the value without a drug metabolite
        FBAorg=optimizeCbModel(model,'max');
        fluxesTmp{i}{1}=FBAorg;
        
        for j = 1:length(drugExchanges)
            model = useDiet(model,basicCompounds);
            model.lb(find(strncmp(model.rxns,'sink_',5)))=0;
            
            if ~isempty(find(ismember(model.rxns,drugExchanges{j,2})))
                % provide the drug
                modelExch=changeRxnBounds(model,drugExchanges{j,2},-1,'l');
                FBA=optimizeCbModel(modelExch,'max');
                fluxesTmp{i}{j+1}=FBA;
            else
                fluxesTmp{i}{j+1}={};
            end
        end
    end
    
    % collect the results
    drugPredictions={};
    drugPredictions{1,2}='None';
    for i = 1:length(orgs)
        drugPredictions{i+1,1}=strrep(orgs{i,1},'.mat','');
        FBA=fluxesTmp{i}{1};
        drugPredictions{i+1,2}=FBA.f;
        for j = 1:length(drugExchanges)
            drugPredictions{1,j+2}=[database.metabolites{find(strcmp(database.metabolites(:,1),drugExchanges{j,1})),2} '_' met];

            if ~isempty(fluxesTmp{i}{j+1})
                FBA=fluxesTmp{i}{j+1};
                drugPredictions{i+1,j+2}=FBA.f;
            else
                drugPredictions{i+1,j+2}=drugPredictions{i+1,2};
            end
        end
    end
    % remove the ones not producing anything
    cnt=1;
    delArray=[];
    for j=2:size(drugPredictions,1)
        if abs(sum(cell2mat(drugPredictions(j,4:end))))<0.0001
            delArray(cnt,1)=j;
            cnt=cnt+1;
        end
    end
    drugPredictions(delArray,:)=[];
    cnt=1;
    delArray=[];
    for j=4:size(drugPredictions,2)
        if abs(sum(cell2mat(drugPredictions(2:end,j))))<0.0001
            delArray(cnt,1)=j;
            cnt=cnt+1;
        end
    end
    drugPredictions(:,delArray)=[];
    drugPredictions=cell2table(drugPredictions);
    writetable(drugPredictions,[resultsFolder filesep 'AGORA2_drugYields_' met],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
end

%% Combine the yields into one file for the supplemental material

fileList = {'AGORA2_drugYields_atp.txt','AGORA2_drugYields_co2.txt','AGORA2_drugYields_nh4.txt','AGORA2_drugYields_pyr.txt'};

orgs={};
for i=1:length(fileList)
    yields = readInputTableForPipeline([resultsFolder fileList{i}]);
    orgs=vertcat(orgs,yields(2:end,1));
end
orgs=unique(orgs);

strainYields={};
for i=1:length(orgs)
    strainYields{i+1,1}=orgs{i};
end

for i=1:length(fileList)
       yields = readInputTableForPipeline([resultsFolder fileList{i}]);
    if size(yields,1) >1
        cols=size(strainYields,2);
        strainYields(1,cols+1:cols+size(yields,2)-1)=yields(1,2:end);
        strainYields(2:end,cols+1:cols+size(yields,2)-1)={'0'};
        for j=2:size(yields,1)
            findStr=find(strcmp(strainYields(:,1),yields{j,1}));
            strainYields(findStr,cols+1:cols+size(yields,2)-1)=yields(j,2:end);
        end
    end
end
writetable(cell2table(strainYields),[resultsFolder filesep 'Table_S8.txt'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
