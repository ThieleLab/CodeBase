
%% Compute conversion capability for Figure 4c

global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER)
    initCobraToolbox
end
solver = CBT_LP_SOLVER;

if numWorkers>0 && ~isempty(ver('parallel'))
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end
environment = getEnvironment();

resultsFolder = [rootDir filesep 'ComputedDrugFluxes'];

database=loadVMHDatabase;

WesternDiet = readtable('WesternDietAGORA2.txt', 'Delimiter', '\t');
WesternDiet=table2cell(WesternDiet);
WesternDiet=cellstr(string(WesternDiet));

taxonomy = readInputTableForPipeline('AGORA2_infoFile.xlsx');

%% compute conversion potential for Figure 4c)

drugExchanges={'InputMet','InputExchange','OutputMet','OutputExchange','Reaction';'glcur','EX_glcur(e)',[],[],'Glucuronic acid';'sn38g','EX_sn38g(e)','sn38','EX_sn38(e)','beta-glucuronidase';'dfdcytd','EX_dfdcytd(e)','dfduri','EX_dfduri(e)','Cytidine deaminase';'fcsn','EX_fcsn(e)','5fura','EX_5fura(e)','Cytosine deaminase';'5fura','EX_5fura(e)','dh5fura','EX_dh5fura(e)','Dihydrouracil dehydrogenase';'4hphac','EX_4hphac(e)','pcresol','EX_pcresol(e)','4-hydroxyphenylacetate decarboxylase';'r788','EX_r788(e)','r406','EX_r406(e)','Alkaline phosphatase';'bzd','EX_bzd(e)','5asa','EX_5asa(e)','Azoreductase';'lactl','EX_lactl(e)','gal','EX_gal(e)','beta-galactosidase';'chlphncl','EX_chlphncl(e)','nchlphncl','EX_nchlphncl(e)','Nitroreductase';'5asa','EX_5asa(e)','ac5asa','EX_ac5asa(e)','Arylamine N-acetyltransferase';'digoxin','EX_digoxin(e)','dihydro_digoxin','EX_dihydro_digoxin(e)','Cardiac glycoside reductase';'srv','EX_srv(e)','bvu','EX_bvu(e)','Pyrimidine-nucleoside phosphorylase';'34dhphe','EX_34dhphe(e)','dopa','EX_dopa(e)','Tryptophan decarboxylase';'tchola','EX_tchola(e)','cholate','EX_cholate(e)','Bile salt hydrolase';'dopa','EX_dopa(e)','mtym','EX_mtym(e)','Dopamine dehydroxylase'};

% fill in the table
drugPredictions={};
drugPredictions{1,1}='ModelID';
for j = 2:length(drugExchanges)
    j
    drugPredictions{1,j}=drugExchanges{j,5};
    fluxesTmp={};
    % parallelisation for faster computation
    parfor i = 2:size(taxonomy,1)
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);

        model=readCbModel([refinedFolder filesep taxonomy{i,1} '.mat']);
        model.lb(find(strncmp(model.rxns,'sink_',5)))=-1;
        model = useDiet(model,WesternDiet);
        model=changeRxnBounds(model,'EX_o2(e)',-10,'l');
        % model = useDiet(model,basicCompounds);
        if ~isempty(find(ismember(model.rxns,drugExchanges{j,2})))
            modelExch=changeRxnBounds(model,drugExchanges{j,2},-1,'l');
            modelExch=changeObjective(modelExch,drugExchanges{j,2});
            FBA=optimizeCbModel(modelExch,'min');
            fluxesTmp{i}=FBA;
        else
            fluxesTmp{i}=[];
        end
    end
    for i = 2:size(taxonomy,1)
        drugPredictions{i,1} = taxonomy{i,1};
        FBA=fluxesTmp{i};
        if ~isempty(FBA)
            if abs(FBA.f) > 0.1
                drugPredictions{i,j}=1;
            else
                drugPredictions{i,j}=0;
            end
        else
            drugPredictions{i,j}=0;
        end
    end
    save([resultsFolder filesep 'drugPredictions.mat'],'drugPredictions');
end

%% summarize drug conversion capabilities by phylum-Figure 4c

load([resultsFolder filesep 'drugPredictions.mat'])

% get organisms with comparative genomics
genomeAnnotation = readtable('gapfilledGenomeAnnotation.txt', 'ReadVariableNames', false, 'Delimiter', 'tab','TreatAsEmpty',['UND. -60001','UND. -2011','UND. -62011']);
genomeAnnotation = table2cell(genomeAnnotation);
orgs=unique(genomeAnnotation(:,1));

drugPredictions(:,2)=[];
[C,IA] = setdiff(drugPredictions(:,1),orgs,'stable');
drugPredictions(IA(2:end),:)=[];

% get taxonomic assignments
taxonomy_reduced=taxonomy;
[C,IA] = setdiff(taxonomy_reduced(:,1),drugPredictions(:,1),'stable');
taxonomy_reduced(IA(2:end),:)=[];
taxonomy_reduced(:,4:10)=strrep(taxonomy_reduced(:,4:10),',','_');
taxCol=find(strcmp(taxonomy_reduced(1,:),'Phylum'));
get_phyla=unique(taxonomy_reduced(2:end,taxCol));
capabilities={};
for j=2:size(drugPredictions,2)
    capabilities{1,j}=drugPredictions{1,j};
end
capabilities{1,1}='Phylum';
for i=1:length(get_phyla)
    capabilities{i+1,1}=get_phyla{i};
    capabilities(i+1,2:end)={0};
    for j=2:size(drugPredictions,2)
        capabilities{i+1,j}=capabilities{i+1,j};
    end
    alltax=taxonomy_reduced(find(strcmp(taxonomy_reduced(:,taxCol),get_phyla{i})),1);
    for k=1:length(alltax)
        findtax=find(strcmp(drugPredictions(:,1),alltax{k}));
        for j=2:size(drugPredictions,2)
            capabilities{i+1,j}=capabilities{i+1,j} + drugPredictions{findtax,j};
        end
    end
end

capabilities(1,:)=strrep(capabilities(1,:),'-','_');
capabilities(1,:)=strrep(capabilities(1,:),' ','_');
capabilities=cell2table(capabilities);
writetable(capabilities,[resultsFolder filesep 'AGORA2_Drug_PhylumSummaries'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
