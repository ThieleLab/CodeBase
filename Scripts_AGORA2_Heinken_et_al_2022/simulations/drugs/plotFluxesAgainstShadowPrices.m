%% plot reaction abundances vs. production potential with shadow prices

% create folder for the plots
spPlotsPath = [rootDir filesep 'Modeling_CRC' filesep 'Plots' filesep 'ShadowPrices'];
mkdir(spPlotsPath)

% load reaction abundances
load([rootDir filesep 'Modeling_CRC' filesep 'Plots' filesep 'ReactionAbundances_Drugs.mat'])

% Japanese diet
fluxPath = [rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices_JD' filesep 'AGORA2_CRC_Objectives_JD.txt'];
% load fluxes
data=load([solutionFolder filesep 'AGORA2_CRC_Objectives_JD.mat']);
fluxes=table2cell(data.('objectives'));
fluxes(:,1)=strrep(fluxes(:,1),'microbiota_model_diet_','');

% Shadow prices
spPath = [rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices_JD' filesep 'AGORA2_CRC_ShadowPrices_JD.txt'];
% load fluxes
shadowPrices = readInputTableForPipeline(spPath);
shadowPrices(1,:)=strrep(shadowPrices(1,:),'microbiota_model_diet_','');

fluxLabels={'EX_sn38[fe]','Diet_EX_sn38g[d]','Deglucuronidated irinotecan from irinotecan';'EX_r406[fe]','Diet_EX_r788[d]','R406 from R788 (fostamatinib)';'EX_5fura[fe]','Diet_EX_fcsn[d]','5-fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_fcsn[d]','56-Dihydro-5-Fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_5fura[d]','56-Dihydro-5-Fluorouracil from 5-fluorouracil';'EX_dfduri[fe]','Diet_EX_dfdcytd[d]','2''2''-Difluorodeoxyuridine from gemcitabine';'EX_dihydro_digoxin[fe]','Diet_EX_digoxin[d]','Dihydrodigoxin from digoxin';'EX_ac5asa[fe]','Diet_EX_5asa[d]','N-acetyl-5-aminosalicylic acid from 5-aminosalicylic acid';'EX_5asa[fe]','Diet_EX_bzd[d]','5-aminosalicylic acid from balsalazide';'EX_ac5asa[fe]','Diet_EX_bzd[d]','N-acetyl-5-aminosalicylic acid from balsalazide';'EX_nchlphncl[fe]','Diet_EX_chlphncl[d]','Nitrosochloramphenicol from chloramphenicol';'EX_bvu[fe]','Diet_EX_srv[d]','(E)-5-(2-Bromovinyl)Uracil from sorivudine';'EX_dopa[fe]','Diet_EX_34dhphe[d]','Dopamine from levodopa';'EX_mtym[fe]','Diet_EX_34dhphe[d]','m-Tyramine from levodopa';'EX_pcresol[fe]','Diet_EX_4hphac[d]','p-cresol from 4-hydroxyphenylacetate';'EX_cholate[fe]','Diet_EX_tchola[d]','Cholic acid from taurocholic acid'};

% pair fluxes and reactions
toPlotWith={'EX_sn38[fe]','SN38G_GLCAASE';'EX_r406[fe]','ALKP_R788';'EX_5fura[fe]','FCSND';'EX_dh5fura[fe]','FURADH';'EX_dfduri[fe]','DFDCYTDD';'EX_dihydro_digoxin[fe]','DIHYDRO_DIGOXINc';'EX_5asa[fe]','BZD_AR_NADi';'EX_ac5asa[fe]','AC5ASA';'EX_nchlphncl[fe]','CZP_NR';'EX_bvu[fe]','PYNP_BRV';'EX_dopa[fe]','3HLYTCL';'EX_mtym[fe]','DOPADH';'EX_pcresol[fe]','4HPHACDC';'EX_cholate[fe]','TCHOLBHS'};

cd(spPlotsPath)

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

cd(rootDir)
