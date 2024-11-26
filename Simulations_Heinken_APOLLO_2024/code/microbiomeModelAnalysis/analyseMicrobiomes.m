
% extract the data for each defined scenario (comparison of microbiome
% samples) and perform simulations

clear all
rootDir = pwd;

defineScenarios

datasets={
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'ModelStatistics.csv'],'Model_sizes'
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'normalizedAbundance.csv'],'Organism_abundance'
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'ReactionAbundance.csv'],'Reaction_abundance'
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'SubsystemAbundance.csv'],'Subsystem_abundance'
    [rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'ReactionPresence.csv'],'Reaction_presence'
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
        writetable(cell2table(data),[datasets{j,2} '_' scenarios{i} '.csv'],'writeVariableNames',false)
    end
    
%     if any(strcmp(scenarios{i},{'IBD_vs_healthy','Infants_undernourished_vs_healthy','PD_vs_healthy'}))
%         % predict a set of metabolites for three scenarios
%         initCobraToolbox
%         solverOK=changeCobraSolver('ibm_cplex','LP');
%         numWorkers=12;
%         computeProfiles=false;
%         
%         resPath = [rootDir filesep 'data' filesep 'GutMicrobiomes' filesep 'MicrobiomeModels'];
%         dietFilePath=[rootDir filesep 'input' filesep 'AverageEuropeanDiet_modified'];
%         abunFilePath = ['Organism_abundance_' scenarios{i} '.csv'];
%         infoFilePath = [scenarios{i} '_samples.csv'];
%         [init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] = initMgPipe(modPath, abunFilePath, computeProfiles, 'resPath', resPath, 'numWorkers', numWorkers, 'dietFilePath', dietFilePath);
%         % compute selected metabolites of interest
%         [objectives,shadowPrices]=analyseObjectiveShadowPrices([resPath filesep 'Diet'], metList, 'numWorkers', numWorkers);
%         writetable(cell2table(objectives),[pwd filesep 'Objectives_AED'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
%     end
    cd ..
end
cd(rootDir)
