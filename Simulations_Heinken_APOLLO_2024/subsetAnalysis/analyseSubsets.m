
datasets={
    'Adults_body_sites_healthy' % nasal cavity, skin, and vagina samples
    'Adults_vs_infants_healthy' % all healthy adults and infants
    'IBD_vs_healthy' % from PMID:24629344
    'Infants_premature_vs_healthy' % all healthy and premature infants
    'Infants_undernourished_vs_healthy' % undernourished and normal infants from Bangladesh
    'Infection_antibiotics_vs_no_antibiotics' % Cholera study, no REF
    'Infection_resistant_vs_susceptible' % from PMID:30057943
    'Infection_vs_healthy' % all healthy gut samples vs. infection
    'Obesity_vs_normalweight' % from PMID:23985870, other samples with BMI available
    'PD_vs_healthy' % from PMID:28662719
    'T2D_vs_healthy' % all T2D vs healthy adults
    };

initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
% 
% for modeling-some potentially interesting metabolites
metList={'EX_ac[fe]'
    'EX_ppa[fe]'
    'EX_but[fe]'
    'EX_for[fe]'
    'EX_isobut[fe]'
    'EX_isoval[fe]'
    'EX_lac_D[fe]'
    'EX_lac_L[fe]'
    'EX_etoh[fe]'
    'EX_succ[fe]'
    'EX_h2[fe]'
    'EX_ch4[fe]'
    'EX_h2s[fe]'
    'EX_tma[fe]'
    'EX_12dhchol[fe]'
    'EX_dchac[fe]'
    'EX_HC02191[fe]'
    'EX_HC02194[fe]'
    'EX_4hphac[fe]'
    'EX_pppn[fe]'
    'EX_ind3ac[fe]'
    'EX_ind3ppa[fe]'
    'EX_leu_L[fe]'
    'EX_phe_L[fe]'
    'EX_tyr_L[fe]'
    'EX_trp_L[fe]'
    'EX_indole[fe]'
    'EX_pcresol[fe]'
    'EX_leu_L[fe]'
    'EX_ile_L[fe]'
    'EX_val_L[fe]'};

numWorkers=20;

combinedDatasets={
    [pwd filesep 'Analysis_microbiome_models' filesep 'Combined_data' filesep 'ModelStatisticsCombined.csv'],'Model sizes','Number'
    [pwd filesep 'Analysis_microbiome_models' filesep 'Combined_data' filesep 'normalizedAbundanceCombined.csv'],'Organism abundance','Relative abundance'
    [pwd filesep 'Analysis_microbiome_models' filesep 'Combined_data' filesep 'ReactionAbundanceCombined.csv'],'Reaction abundance','Relative abundance'
    [pwd filesep 'Analysis_microbiome_models' filesep 'Combined_data' filesep 'SubsystemAbundanceCombined.csv'],'Subsystem abundance','Relative abundance'
    [pwd filesep 'Analysis_microbiome_models' filesep 'Combined_data' filesep 'ReactionPresenceCombined.csv'],'Reactions presence','Presence'
    };

cd([pwd filesep 'Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups'])

for i=1:length(datasets)
    cd(datasets{i})

    metadata=readInputTableForPipeline([datasets{i} '_samples.csv']);

    for j=1:size(combinedDatasets,1)
        j
        dataTable = readInputTableForPipeline(combinedDatasets{j,1});
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
        cell2csv([strrep(combinedDatasets{j,2}, ' ','_') '_' datasets{i} '.csv'],data)

        if strcmp(datasets{i},'Infection_antibiotics_vs_no_antibiotics')
            strat='Antibiotics';
        elseif strcmp(datasets{i},'Infection_resistant_vs_susceptible')
            strat='Stratification';
        elseif strcmp(datasets{i},'Adults_body_sites_healthy')
            strat='Body site';
        elseif strcmp(datasets{i},'Adults_vs_infants_healthy')
            strat='Age group';
        else
            strat='Disease name';
        end

        if j==1 || j==4
            mkdir('ViolinPlots')
            cd('ViolinPlots')

            makeViolinPlots(data, metadata, 'stratification', strat, 'plottedFeature',combinedDatasets{j,2},'unit',combinedDatasets{j,3});
            cd ..
        end

%         if j==3
%             mkdir('ViolinPlots_Reactions')
%             cd('ViolinPlots_Reactions')
% 
%             makeViolinPlots(data, metadata, 'stratification', strat, 'plottedFeature',combinedDatasets{j,2},'unit',combinedDatasets{j,3});
%             cd ..
%         end

        if j==1 || j==5
            mkdir('StatisticalAnalysis')
            cd('StatisticalAnalysis')
            [Statistics,significantFeatures] = performStatisticalAnalysis(data', metadata, 'stratification', strat);
            Statistics(:,2)=strrep(Statistics(:,2),',',' ');
            cell2csv(['Statistics_' strrep(combinedDatasets{j,2}, ' ','_') '.csv'],Statistics)
            cell2csv(['SignificantFeatures_' strrep(combinedDatasets{j,2}, ' ','_') '.csv'],significantFeatures)
            cd ..
        end
    end

%     perform microbiome modeling-only for some subsets, would take very long otherwise
    if i==3 || i==10
        % path to microbiome models
        origModPath='C:\Users\Almut\OneDrive - National University of Ireland, Galway\90000_genomes_Almeida_2019\MicrobiomeModelling\Microbiome_models_AED';

        mkdir('Models')
        % copy all of the respective microbiome models
        for k=2:size(metadata,1)
            copyfile([origModPath filesep 'microbiota_model_samp_' metadata{k,1} '.mat'],['Models' filesep 'microbiota_model_samp_' metadata{k,1} '.mat'])
        end
        modPath = 'Models';

        mkdir(['MicrobiomeModeling_' datasets{i}])
        spPath = ['MicrobiomeModeling_' datasets{i}];

        [objectives,shadowPrices]=analyseObjectiveShadowPrices(modPath, metList, 'resultsFolder', spPath, 'numWorkers', numWorkers);

        objectives=cell2table(objectives);
        writetable(objectives,'Objectives_AED','FileType','text','WriteVariableNames',false,'Delimiter','tab');
    end
    cd ..
end
