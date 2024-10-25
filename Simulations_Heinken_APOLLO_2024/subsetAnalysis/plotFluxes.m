
metNames = {'EX_ac[fe]','Acetate';'EX_ppa[fe]','Propionate';'EX_but[fe]','Butyrate';'EX_isobut[fe]','Isobutyrate';'EX_isoval[fe]','Isovalerate';'EX_lac_D[fe]','D-lactate';'EX_lac_L[fe]','L-lactate';'EX_for[fe]','Formate';'EX_etoh[fe]','Ethanol';'EX_h2s[fe]','Hydrogen sulfide';'EX_tma[fe]','Trimethylamine';'EX_phe_L[fe]','Phenylalanine';'EX_tyr_L[fe]','Tyrosine';'EX_trp_L[fe]','Tryptophan';'EX_dopa[fe]','Dopamine';'EX_taur[fe]','Taurine';'EX_hcys_L[fe]','Homocysteine';'EX_pcresol[fe]','p-cresol';'EX_ind3ac[fe]','Indoleacetate';'EX_indole[fe]','Indole';'EX_4abut[fe]','GABA';'EX_leu_L[fe]','Leucine';'EX_ile_L[fe]','Isoleucine';'EX_val_L[fe]','Valine';'EX_cholate[fe]','Cholate';'EX_dchac[fe]','Deoxycholate';'EX_HC02191[fe]','Lithocholate';'EX_HC02194[fe]','Ursodeoxycholate';'EX_12dhchol[fe]','12-dehydrocholate';'EX_7ocholate[fe]','7-dehydrocholate'};

mkdir(['Analysis_microbiome_models' filesep 'Summary_for_figures' filesep 'FluxPlots'])


datasets={
    'IBD_vs_healthy' % from PMID:24629344
    'Infants_undernourished_vs_healthy' % undernourished and normal infants from Bangladesh
    'PD_vs_healthy' % from PMID:28662719
    };

% create boxplots for fluxes where it applies

for d=1:length(datasets)
    fluxes = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{d} filesep 'Objectives.txt']);
    fluxes(1,:) = strrep(fluxes(1,:),'microbiota_model_samp_','');
    % remove nonsignificant fluxes
    statResults = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{d} filesep datasets{d} '_stat_flux.csv']);
    statResults(1,:) = [];
    findNS = find(cell2mat(statResults(:,end))>=0.05);
    statResults(findNS,:) = [];
    findNS = find(isnan(cell2mat(statResults(:,end))));
    statResults(findNS,:) = [];
    [C,I] = setdiff(fluxes(:,1),statResults(:,1),'stable');
    fluxes(I(2:end),:) = [];

    mets={};
    for i=2:size(fluxes,1)
        mets{i} = metNames{find(strcmp(metNames(:,1),fluxes{i,1})),2};
    end

        % cut to 12 metabolites if there are more
        if size(fluxes,1)>12
            summed = sum(cell2mat(fluxes(2:end,2:end)),2);
            summed = vertcat(1000000,summed);
            [I,J] = sort(summed,'descend');
            fluxes(J(14:end),:) = [];
        end

    samples = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Subgroup_analysis' filesep 'Subgroups' filesep datasets{d} filesep datasets{d} '_samples.csv']);
    stratCol = find(strcmp(samples(1,:),'Disease name'));
    samples(:,stratCol) = strrep(samples(:,stratCol),'Parkinson','PD');
    if d==1
        stratGroups = {'Healthy','CD','UC'};
        cols = [1 1 0
            1 0 0
            0 0 1];
    elseif d==2
        stratGroups = {'Healthy','Undernourished'};
        cols = [1 0 0
            0 0 1];
    elseif d==3
        stratGroups = {'Healthy','PD'};
        cols = [1 0 0
            0 0 1];
    end

    Groups={};
    for i=2:size(fluxes,2)
        Groups{i-1,1}=samples{find(strcmp(samples(:,1),fluxes{1,i})),stratCol};
    end

    %     table = {};
    %     cnt = 1;
    %     for i=2:size(fluxes,1)
    %         for j=2:size(fluxes,2)
    %             table{cnt,1} = fluxes{i,j};
    %             table{cnt,2} = samples(find(strcmp(samples(:,1),fluxes{1,j})),stratCol);
    %             table{cnt,3} = mets{i};
    %             cnt=cnt+1;
    %         end
    %     end

    %     table=cell2table(table,'VariableNames',{'Value','Groups','Metabolite'});
    %     table.Metabolite=categorical(table.Metabolite,mets(2:end));
    Groups=categorical(Groups,stratGroups);
    figure
    if d==1
        for i=2:size(fluxes,1)
            %     boxchart(table.Metabolite,table.Value,'GroupByColor',table.Groups)
            subplot(3,4,i-1)
            boxchart(cell2mat(fluxes(i,2:end)),'GroupByColor',Groups)
            legend(stratGroups,'Location','best')
            xticklabels(mets{i})
            set(gca, 'FontSize',14)
            ylabel('Flux (mmol/person/day)')
        end
    elseif d==2
        for i=2:size(fluxes,1)
            %     boxchart(table.Metabolite,table.Value,'GroupByColor',table.Groups)
            subplot(3,3,i-1)
            boxchart(cell2mat(fluxes(i,2:end)),'GroupByColor',Groups)
            legend(stratGroups,'Location','best')
            xticklabels(mets{i})
            set(gca, 'FontSize',14)
            ylabel('Flux (mmol/person/day)')
        end
    elseif d==3
        for i=2:size(fluxes,1)
            %     boxchart(table.Metabolite,table.Value,'GroupByColor',table.Groups)
            subplot(1,3,i-1)
            boxchart(cell2mat(fluxes(i,2:end)),'GroupByColor',Groups)
            legend(stratGroups,'Location','best')
            xticklabels(mets{i})
            set(gca, 'FontSize',14)
            ylabel('Flux (mmol/person/day)')
        end
    end
end
