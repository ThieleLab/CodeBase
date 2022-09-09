%% plot reaction abundances vs. production potential

solutionFolder = [rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices_JD'];
plotsPath = [rootDir filesep 'Modeling_CRC' filesep 'Plots'];
mkdir(plotsPath)

% load the predicted metabolite fluxes
abundancePath=[rootDir filesep 'normalizedCoverage.csv'];
abundance = readInputTableForPipeline(abundancePath);
modelPath=[rootDir filesep 'Modeling_CRC' filesep 'Species'];

rxnLabels={'DFDCYTDD','Cytidine deaminase';'FCSND','Cytosine deaminase';'FURADH','Dihydropyrimidine dehydrogenase';'ALKP_R788','Alkaline phosphatase';'BZD_AR_NADi','Azoreductase';'CZP_NR','Nitroreductase';'AC5ASA','Arylamine N-acetyltransferase';'DIHYDRO_DIGOXINc','Digoxin reductase';'SN38G_GLCAASE','Beta-glucuronidase';'PYNP_BRV','Pyrimidine-nucleoside phosphorylase';'3HLYTCL','3-Hydroxy-L-Tyrosine Carboxy-Lyase';'DOPADH','Dopamine dehydroxylase';'4HPHACDC','4-hydroxyphenylacetate decarboxylase';'TCHOLBHS','Bile salt hydrolase'};

% summarize names to get all individual reactions
rxnsToSummarize={'FCSNDepp','FCSND';'FCSNDe','FCSND';'FCSND','FCSND';'CZP_NRepp','CZP_NR';'CZP_NRe','CZP_NR';'CZP_NR','CZP_NR';'AC5ASAc','AC5ASA';'AC5ASAepp','AC5ASA';'AC5ASAe','AC5ASA';'SN38G_GLCAASEepp','SN38G_GLCAASE';'SN38G_GLCAASEe','SN38G_GLCAASE';'SN38G_GLCAASE','SN38G_GLCAASE';'DFDCYTDDepp','DFDCYTDD';'DFDCYTDDe','DFDCYTDD';'DFDCYTDD','DFDCYTDD';'ALKP_R788epp','ALKP_R788';'ALKP_R788e','ALKP_R788';'ALKP_R788','ALKP_R788';'56DFURADHepp','FURADH';'56DFURADHe','FURADH';'5FURADH','FURADH';'56DFURADH','FURADH';'TCHOLBHSepp','TCHOLBHS';'TCHOLBHSe','TCHOLBHS';'TCHOLAH','TCHOLBHS';'DIHYDRO_DIGOXINc','DIHYDRO_DIGOXINc';'PYNP_BRV','PYNP_BRV';'3HLYTCL','3HLYTCL';'DOPADH','DOPADH';'4HPHACDC','4HPHACDC';'BZD_AR_NADi','BZD_AR_NADi'};

reactions={};
% recalculate abundances on a species to species basis
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
save([rootDir filesep 'Modeling_CRC' filesep 'Plots' filesep 'ReactionAbundances_Drugs.mat'],'reactions')

% load fluxes
data=load([solutionFolder filesep 'AGORA2_CRC_Objectives_JD.mat']);
fluxes=table2cell(data.('objectives'));
fluxes(:,1)=strrep(fluxes(:,1),'microbiota_model_diet_','');

fluxLabels={'EX_sn38[fe]','Diet_EX_sn38g[d]','Deglucuronidated irinotecan from irinotecan';'EX_r406[fe]','Diet_EX_r788[d]','R406 from R788 (fostamatinib)';'EX_5fura[fe]','Diet_EX_fcsn[d]','5-fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_fcsn[d]','56-Dihydro-5-Fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_5fura[d]','56-Dihydro-5-Fluorouracil from 5-fluorouracil';'EX_dfduri[fe]','Diet_EX_dfdcytd[d]','2''2''-Difluorodeoxyuridine from gemcitabine';'EX_dihydro_digoxin[fe]','Diet_EX_digoxin[d]','Dihydrodigoxin from digoxin';'EX_ac5asa[fe]','Diet_EX_5asa[d]','N-acetyl-5-ASA from 5-ASA';'EX_5asa[fe]','Diet_EX_bzd[d]','5-ASA from balsalazide';'EX_ac5asa[fe]','Diet_EX_bzd[d]','N-acetyl-5-ASA from balsalazide';'EX_nchlphncl[fe]','Diet_EX_chlphncl[d]','Nitrosochloramphenicol from chloramphenicol';'EX_bvu[fe]','Diet_EX_srv[d]','(E)-5-(2-Bromovinyl)Uracil from sorivudine';'EX_dopa[fe]','Diet_EX_34dhphe[d]','Dopamine from levodopa';'EX_mtym[fe]','Diet_EX_34dhphe[d]','m-Tyramine from levodopa';'EX_pcresol[fe]','Diet_EX_4hphac[d]','p-cresol from 4-hydroxyphenylacetate';'EX_cholate[fe]','Diet_EX_tchola[d]','Cholic acid from taurocholic acid'};

% load abundances
load([rootDir filesep 'Modeling_CRC' filesep 'Plots' filesep 'ReactionAbundances_Drugs.mat'])

% pair fluxes and reactions
toPlotWith={'EX_sn38[fe]','SN38G_GLCAASE';'EX_r406[fe]','ALKP_R788';'EX_5fura[fe]','FCSND';'EX_dh5fura[fe]','FURADH';'EX_dfduri[fe]','DFDCYTDD';'EX_dihydro_digoxin[fe]','DIHYDRO_DIGOXINc';'EX_5asa[fe]','BZD_AR_NADi';'EX_ac5asa[fe]','AC5ASA';'EX_nchlphncl[fe]','CZP_NR';'EX_bvu[fe]','PYNP_BRV';'EX_dopa[fe]','3HLYTCL';'EX_mtym[fe]','DOPADH';'EX_pcresol[fe]','4HPHACDC';'EX_cholate[fe]','TCHOLBHS'};

%% plot fluxes against reaction abundances
f=figure;
% the colors are created during the script plotDrugViolins!
for i=2:size(fluxes,2)
    subplot(4,4,i-1)
    data=cell2mat(fluxes(3:end,i));
    % find the matching reaction
    findFlux=find(strcmp(toPlotWith(:,1),fluxes{1,i}));
    findRxn=find(strcmp(reactions(1,:),toPlotWith{findFlux,2}));
    for k=3:size(fluxes,1)
        findInd=find(strcmp(reactions(:,1),fluxes{k,1}));
        data(k-2,2)=str2double(reactions{findInd,findRxn});
    end
    scatter(data(:,1),data(:,2),[],cols(i-1,:),'filled','o','MarkerEdgeColor','black')
    find1=find(strcmp(fluxLabels(:,1),fluxes{1,i}));
    find2=find(strcmp(fluxLabels(:,2),fluxes{2,i}));
    h=xlabel([fluxLabels{intersect(find1,find2),3} ' flux']);
    set(h,'interpreter','none')
    ylabel('Abundance')
    set(gca, 'FontSize', 11);
    h=title(rxnLabels{find(strcmp(rxnLabels(:,1),toPlotWith{findFlux,2})),2});
    set(h,'interpreter','none')
end
f.Renderer='painters';

