%% plot reaction abundances vs. production potential with shadow prices

% load the predicted metabolite fluxes
plotsPath = [rootDir filesep 'Modeling_CRC' filesep 'Plots' filesep 'ShadowPrices'];
mkdir(plotsPath)

abundancePath=[rootDir filesep 'Modeling_CRC' filesep 'normCoverage_CRC.csv'];
abundance = readInputTableForPipeline(abundancePath);
modelPath=[rootDir filesep 'panModelsAGORA2' filesep 'Species'];

rxnLabels={'DFDCYTDD','Cytidine deaminase';'FCSND','Cytosine deaminase';'FURADH','Dihydropyrimidine dehydrogenase';'ALKP_R788','Alkaline phosphatase';'BZD_AR_NAD','Azoreductase';'CZP_NR','Nitroreductase';'AC5ASA','Arylamine N-acetyltransferase';'DIHYDRO_DIGOXINc','Digoxin reductase';'SN38G_GLCAASE','Beta-glucuronidase';'PYNP_BRV','Pyrimidine-nucleoside phosphorylase';'3HLYTCL','3-Hydroxy-L-Tyrosine Carboxy-Lyase';'DOPADH','Dopamine dehydroxylase';'4HPHACDC','4-hydroxyphenylacetate decarboxylase';'TCHOLBHS','Bile salt hydrolase'};

reactions={};
% recalculate abundances on a species to species basis

% summarize names to get all individual reactions
rxnsToSummarize={'FCSNDepp','FCSND';'FCSNDe','FCSND';'FCSND','FCSND';'CZP_NRepp','CZP_NR';'CZP_NRe','CZP_NR';'CZP_NR','CZP_NR';'AC5ASAc','AC5ASA';'AC5ASAepp','AC5ASA';'AC5ASAe','AC5ASA';'SN38G_GLCAASEepp','SN38G_GLCAASE';'SN38G_GLCAASEe','SN38G_GLCAASE';'SN38G_GLCAASE','SN38G_GLCAASE';'DFDCYTDDepp','DFDCYTDD';'DFDCYTDDe','DFDCYTDD';'DFDCYTDD','DFDCYTDD';'ALKP_R788epp','ALKP_R788';'ALKP_R788e','ALKP_R788';'ALKP_R788','ALKP_R788';'56DFURADHepp','FURADH';'56DFURADHe','FURADH';'5FURADH','FURADH';'56DFURADH','FURADH';'TCHOLBHSepp','TCHOLBHS';'TCHOLBHSe','TCHOLBHS';'TCHOLAH','TCHOLBHS';'DIHYDRO_DIGOXINc','DIHYDRO_DIGOXINc';'PYNP_BRV','PYNP_BRV';'3HLYTCL','3HLYTCL';'DOPADH','DOPADH';'4HPHACDC','4HPHACDC';'BZD_AR_NAD','BZD_AR_NAD'};

uniqueRxns=unique(rxnsToSummarize(:,2));
for j=1:length(uniqueRxns)
    reactions{1,j+1}=uniqueRxns{j};
end

for i=2:size(abundance,2)
    i
    reactions{i,1}=abundance{1,i};
    for j=1:length(uniqueRxns)
        reactions{i,j+1}={'0'};
    end
    for k=2:size(abundance,1)
        if abundance{k,i} > 0.0000000001
            % load species that are in the sample
            load([modelPath filesep abundance{k,1}]);
            
            % loop through all reaction to analyse
            for j=1:length(uniqueRxns)
                getRxns=rxnsToSummarize(find(strcmp(rxnsToSummarize(:,2),uniqueRxns{j})),1);
                % find if the reactions is in the species
                if ~isempty(intersect(model.rxns,getRxns))
                    reactions{i,j+1}=num2str(str2double(reactions{i,j+1})+abundance{k,i});
                end
            end
        end
    end
end

% Japanese diet
fluxPath = [rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices_JD' filesep 'AGORA2_CRC_Objectives_JD.txt'];
% load fluxes
fluxes = readInputTableForPipeline(fluxPath);
fluxes(:,1)=strrep(fluxes(:,1),'microbiota_model_diet_','');

% Shadow prices
spPath = [rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices_JD' filesep 'AGORA2_CRC_ShadowPrices_JD.txt'];
% load fluxes
shadowPrices = readInputTableForPipeline(spPath);
shadowPrices(1,:)=strrep(shadowPrices(1,:),'microbiota_model_diet_','');

fluxLabels={'EX_sn38[fe]','Diet_EX_sn38g[d]','Deglucuronidated irinotecan from irinotecan';'EX_r406[fe]','Diet_EX_r788[d]','R406 from R788 (fostamatinib)';'EX_5fura[fe]','Diet_EX_fcsn[d]','5-fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_fcsn[d]','56-Dihydro-5-Fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_5fura[d]','56-Dihydro-5-Fluorouracil from 5-fluorouracil';'EX_dfduri[fe]','Diet_EX_dfdcytd[d]','2''2''-Difluorodeoxyuridine from gemcitabine';'EX_dihydro_digoxin[fe]','Diet_EX_digoxin[d]','Dihydrodigoxin from digoxin';'EX_ac5asa[fe]','Diet_EX_5asa[d]','N-acetyl-5-aminosalicylic acid from 5-aminosalicylic acid';'EX_5asa[fe]','Diet_EX_bzd[d]','5-aminosalicylic acid from balsalazide';'EX_ac5asa[fe]','Diet_EX_bzd[d]','N-acetyl-5-aminosalicylic acid from balsalazide';'EX_nchlphncl[fe]','Diet_EX_chlphncl[d]','Nitrosochloramphenicol from chloramphenicol';'EX_bvu[fe]','Diet_EX_srv[d]','(E)-5-(2-Bromovinyl)Uracil from sorivudine';'EX_dopa[fe]','Diet_EX_34dhphe[d]','Dopamine from levodopa';'EX_mtym[fe]','Diet_EX_34dhphe[d]','m-Tyramine from levodopa';'EX_pcresol[fe]','Diet_EX_4hphac[d]','p-cresol from 4-hydroxyphenylacetate';'EX_cholate[fe]','Diet_EX_tchola[d]','Cholic acid from taurocholic acid'};

% pair fluxes and reactions
toPlotWith={'EX_sn38[fe]','SN38G_GLCAASE';'EX_r406[fe]','ALKP_R788';'EX_5fura[fe]','FCSND';'EX_dh5fura[fe]','FURADH';'EX_dfduri[fe]','DFDCYTDD';'EX_dihydro_digoxin[fe]','DIHYDRO_DIGOXINc';'EX_5asa[fe]','BZD_AR_NAD';'EX_ac5asa[fe]','AC5ASA';'EX_nchlphncl[fe]','CZP_NR';'EX_bvu[fe]','PYNP_BRV';'EX_dopa[fe]','3HLYTCL';'EX_mtym[fe]','DOPADH';'EX_pcresol[fe]','4HPHACDC';'EX_cholate[fe]','TCHOLBHS'};

cd(plotsPath)

%% plot zero vs. nonzero shadow prices
for i=2:size(fluxes,2)
    data=cell2mat(fluxes(3:end,i));
    % find the matching reaction
    findFlux=find(strcmp(toPlotWith(:,1),fluxes{1,i}));
    findRxn=find(strncmp(reactions(1,:),toPlotWith{findFlux,2},length(toPlotWith{findFlux,2})));
    for k=3:size(fluxes,1)
        data(k-2,2)=0;
        for j=1:length(findRxn)
            data(k-2,2) = data(k-2,2)+str2double(reactions{find(strcmp(reactions(:,1),fluxes{k,1})),findRxn(j)});
        end
    end
    find1=find(strcmp(shadowPrices(:,2),fluxes{1,i}));
    find2=find(strcmp(shadowPrices(:,3),fluxes{2,i}));
    sps=intersect(find1,find2);
    for j=1:length(sps)
        dataSP=[];
        dataNoSP=[];
        for k=3:size(fluxes,1)
            findInd=find(strcmp(shadowPrices(1,:),fluxes{k,1}));
            if abs(cell2mat(shadowPrices(sps(j),findInd))) > 100
                dataSP(size(dataSP,1)+1,:)=data(k-2,:);
            else
                dataNoSP(size(dataNoSP,1)+1,:)=data(k-2,:);
            end
        end
        if size(dataSP,1) > 4
            f=figure;
            hold on
            if ~isempty(dataNoSP)
                scatter(dataNoSP(:,1),dataNoSP(:,2),'b','filled','o','MarkerEdgeColor','black')
            end
            hold on
            scatter(dataSP(:,1),dataSP(:,2),'m','filled','o','MarkerEdgeColor','black')
            find1=find(strcmp(fluxLabels(:,1),fluxes{1,i}));
            find2=find(strcmp(fluxLabels(:,2),fluxes{2,i}));
            h=xlabel(fluxLabels{intersect(find1,find2),3});
            set(h,'interpreter','none')
            h=ylabel(rxnLabels{find(strcmp(rxnLabels(:,1),toPlotWith{findFlux,2})),2});
            set(h,'interpreter','none')
            set(gca, 'FontSize', 10);
            if ~isempty(dataNoSP)
            xlim([0 (max(vertcat(dataSP(:,1),dataNoSP(:,1))) + max(vertcat(dataSP(:,1),dataNoSP(:,1)))/10)]);
            ylim([0 (max(vertcat(dataSP(:,2),dataNoSP(:,2))) + max(vertcat(dataSP(:,2),dataNoSP(:,2)))/10)]);
            legend('Shadow price zero','Shadow price nonzero','Location','Northwest');
            else
                legend('Shadow price nonzero','Location','Northwest');
            end
            f.Renderer='painters';
            h=title(shadowPrices{sps(j),1});
            set(h,'interpreter','none')
            plotname=strrep([reactions{1,find(strcmp(reactions(1,:),toPlotWith{findFlux,2}))} '_',fluxes{1,i} '_',fluxes{2,i}  '_' shadowPrices{sps(j),1}],'[','');
            plotname=strrep(plotname,']','');
            print(plotname,'-dpng','-r300')
        end
    end
    close all
end
