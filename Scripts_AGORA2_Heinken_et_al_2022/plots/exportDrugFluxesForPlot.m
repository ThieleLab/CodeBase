
% Plot correlations between reaction abundances and fluxes
reactionPath = [rootDir filesep 'Modeling_CRC' filesep 'MicrobiomeModels' filesep 'ReactionAbundance.csv'];
plotsPath = [rootDir filesep 'Modeling_CRC' filesep 'Plots'];
mkdir(plotsPath)

rxnLabels={'DFDCYTDD','Cytidine deaminase';'FCSND','Cytosine deaminase';'FURADH','Dihydropyrimidine dehydrogenase';'ALKP_R788','Alkaline phosphatase';'BZD_AR_NAD','Azoreductase';'CZP_NR','Nitroreductase';'AC5ASA','Arylamine N-acetyltransferase';'DIHYDRO_DIGOXINc','Digoxin reductase';'SN38G_GLCAASE','Beta-glucuronidase';'PYNP_BRV','Pyrimidine-nucleoside phosphorylase';'3HLYTCL','3-Hydroxy-L-Tyrosine Carboxy-Lyase';'DOPADH','Dopamine dehydroxylase';'4HPHACDC','4-hydroxyphenylacetate decarboxylase';'TCHOLBHS','Bile saly hydrolase'};

% load reaction abundances
reactions = readInputTableForPipeline(reactionPath);
% Do not include biomass and sink/demand-not very informative
reactions(find(strncmp(reactions(:,1),'bio',3)),:)=[];
reactions(find(strncmp(reactions(:,1),'DM_',3)),:)=[];
reactions(find(strncmp(reactions(:,1),'sink_',5)),:)=[];

% Japanese diet
fluxPath = [rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices_JD' filesep 'AGORA2_CRC_Objectives_JD.txt'];
% load fluxes
fluxes = readInputTableForPipeline(fluxPath);

fluxLabels={'EX_sn38[fe]','Diet_EX_sn38g[d]','Deglucuronidated irinotecan from irinotecan';'EX_r406[fe]','Diet_EX_r788[d]','R406 from R788 (fostamatinib)';'EX_5fura[fe]','Diet_EX_fcsn[d]','5-fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_fcsn[d]','56-Dihydro-5-Fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_5fura[d]','56-Dihydro-5-Fluorouracil from 5-fluorouracil';'EX_dfduri[fe]','Diet_EX_dfdcytd[d]','2''2''-Difluorodeoxyuridine from gemcitabine';'EX_dihydro_digoxin[fe]','Diet_EX_digoxin[d]','Dihydrodigoxin from digoxin';'EX_ac5asa[fe]','Diet_EX_5asa[d]','N-acetyl-5-aminosalicylic acid from 5-aminosalicylic acid';'EX_5asa[fe]','Diet_EX_bzd[d]','5-aminosalicylic acid from balsalazide';'EX_ac5asa[fe]','Diet_EX_bzd[d]','N-acetyl-5-aminosalicylic acid from balsalazide';'EX_nchlphncl[fe]','Diet_EX_chlphncl[d]','Nitrosochloramphenicol from chloramphenicol';'EX_bvu[fe]','Diet_EX_srv[d]','(E)-5-(2-Bromovinyl)Uracil from sorivudine';'EX_dopa[fe]','Diet_EX_34dhphe[d]','Dopamine from levodopa';'EX_mtym[fe]','Diet_EX_34dhphe[d]','m-Tyramine from levodopa';'EX_pcresol[fe]','Diet_EX_4hphac[d]','p-cresol from 4-hydroxyphenylacetate';'EX_cholate[fe]','Diet_EX_tchola[d]','Cholic acid from taurocholic acid'};

% extract absolute capabilities for heatmap
for i=2:size(fluxes,2)
    find1=find(strcmp(fluxLabels(:,1),fluxes{1,i}));
    find2=find(strcmp(fluxLabels(:,2),fluxes{2,i}));
    fluxes{1,i}=fluxLabels{intersect(find1,find2),3};
end
fluxes(2,:)=[];
fluxes(:,1)=strrep(fluxes(:,1),'microbiota_model_diet_','');
% export quantitative values
cell2csv([plotsPath filesep 'DrugFluxes_CRC_Microbiomes.csv'],fluxes);

% export absolute values
for i=2:size(fluxes,1)
    for j=2:size(fluxes,2)
        if str2double(fluxes{i,j})>0.000001
            fluxes{i,j}=1;
        else
            fluxes{i,j}=0;
        end
    end
end
cell2csv([plotsPath filesep 'absoluteCapabilities.csv'],fluxes);
