
defineScenarios

datasets={
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'ModelStatistics.csv'],'Model sizes','Number'
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'normalizedAbundance.csv'],'Organism abundance','Relative abundance'
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'ReactionAbundance.csv'],'Reaction abundance','Relative abundance'
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'SubsystemAbundance.csv'],'Subsystem abundance','Relative abundance'
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'ReactionPresence.csv'],'Reactions presence','Presence'
    };

% define a set of potentially interesting metabolites
metList={'EX_ac[fe]'
'EX_ppa[fe]'
'EX_but[fe]'
'EX_isobut[fe]'
'EX_isoval[fe]'
'EX_lac_D[fe]'
'EX_lac_L[fe]'
'EX_for[fe]'
'EX_etoh[fe]'
'EX_h2s[fe]'
'EX_tma[fe]'
'EX_phe_L[fe]'
'EX_tyr_L[fe]'
'EX_trp_L[fe]'
'EX_dopa[fe]'
'EX_taur[fe]'
'EX_pcresol[fe]'
'EX_indole[fe]'
'EX_4abut[fe]'
'EX_leu_L[fe]'
'EX_ile_L[fe]'
'EX_val_L[fe]'
'EX_cholate[fe]'
'EX_dchac[fe]'
'EX_HC02191[fe]'
'EX_12dhchol[fe]'
'EX_7ocholate[fe]'};

createScenarios

cd([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Scenarios'])

for i=1:length(scenarios)
    cd(scenarios{i})

    metadata=readInputTableForPipeline([scenarios{i} '_samples.csv']);

    for j=1:size(datasets,1)
        j
        dataTable = readInputTableForPipeline(datasets{j,1});
        if j==2
            data=dataTable;
        else
            data=dataTable';
        end

        % remove samples not in subset
        [C,IA]=setdiff(data(1,:),metadata(2:end,1),'stable');
        data(:,IA(2:end))=[];

        if j==2
            % remove taxa that are not present in the subset
            cnt=1;
            delArray=[];
            for k=2:size(data,1)
                if sum(str2double(data(k,2:end))) < 0.000001
                    delArray(cnt)=k;
                    cnt=cnt+1;
                end
            end
            data(delArray,:) = [];
        end

        % save the subset
        cell2csv([strrep(datasets{j,2}, ' ','_') '_' scenarios{i} '.csv'],data)

        if strcmp(scenarios{i},'Infection_antibiotics_vs_no_antibiotics')
            strat='Antibiotics';
        elseif strcmp(scenarios{i},'Infection_resistant_vs_susceptible')
            strat='Stratification';
        elseif strcmp(scenarios{i},'Adults_body_sites_healthy')
            strat='Body site';
        elseif strcmp(scenarios{i},'Adults_vs_infants_healthy')
            strat='Age group';
        else
            strat='Disease name';
        end

        if j==1 || j==4
            mkdir('ViolinPlots')
            cd('ViolinPlots')

            makeViolinPlots(data, metadata, 'stratification', strat, 'plottedFeature',datasets{j,2},'unit',datasets{j,3});
            cd ..
        end

        if j==1 || j==5
            mkdir('StatisticalAnalysis')
            cd('StatisticalAnalysis')
            [Statistics,significantFeatures] = performStatisticalAnalysis(data', metadata, 'stratification', strat);
            Statistics(:,2)=strrep(Statistics(:,2),',',' ');
            cell2csv(['Statistics_' strrep(datasets{j,2}, ' ','_') '.csv'],Statistics)
            cell2csv(['SignificantFeatures_' strrep(datasets{j,2}, ' ','_') '.csv'],significantFeatures)
            cd ..
        end
    end
    
    if i==3 || i==5 || i==10
        % predict a set of metabolites for three scenarios
        initCobraToolbox
        solverOK=changeCobraSolver('ibm_cplex','LP');
        numWorkers=12;
        computeProfiles=false;
        
        resPath = [rootDir filesep 'data' filesep 'GutMicrobiomes' filesep 'MicrobiomeModels'];
        dietFilePath=[rootDir filesep 'input' filesep 'AverageEuropeanDiet_modified'];
        abunFilePath = ['Organism_abundance_' scenarios{i} '.csv'];
        infoFilePath = [scenarios{i} '_samples.csv'];
        [init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] = initMgPipe(modPath, abunFilePath, computeProfiles, 'resPath', resPath, 'numWorkers', numWorkers, 'dietFilePath', dietFilePath);
        % compute selected metabolites of interest
        [objectives,shadowPrices]=analyseObjectiveShadowPrices([resPath filesep 'Diet'], metList, 'numWorkers', numWorkers);
        writetable(cell2table(objectives),[pwd filesep 'Objectives_AED'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
    end
    cd ..
end
cd(rootDir)
