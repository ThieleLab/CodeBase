
%% plot drug metabolism potential in each CRC model as violin plots
plotsPath = [rootDir filesep 'Modeling_CRC' filesep 'Plots'];

% Japanese diet
fluxPath = [rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices_JD' filesep 'AGORA2_CRC_Objectives_JD.txt'];
% load fluxes
fluxes = readtable(fluxPath, 'ReadVariableNames', false);
fluxes = table2cell(fluxes);

fluxLabels={'EX_sn38[fe]','Diet_EX_sn38g[d]','Deglucuronidated irinotecan from irinotecan';'EX_r406[fe]','Diet_EX_r788[d]','R406 from R788 (fostamatinib)';'EX_5fura[fe]','Diet_EX_fcsn[d]','5-fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_fcsn[d]','56-Dihydro-5-Fluorouracil from 5-fluorocytosine';'EX_dh5fura[fe]','Diet_EX_5fura[d]','56-Dihydro-5-Fluorouracil from 5-fluorouracil';'EX_dfduri[fe]','Diet_EX_dfdcytd[d]','2''2''-Difluorodeoxyuridine from gemcitabine';'EX_dihydro_digoxin[fe]','Diet_EX_digoxin[d]','Dihydrodigoxin from digoxin';'EX_ac5asa[fe]','Diet_EX_5asa[d]','N-acetyl-5-ASA from 5-ASA';'EX_5asa[fe]','Diet_EX_bzd[d]','5-ASA from balsalazide';'EX_ac5asa[fe]','Diet_EX_bzd[d]','N-acetyl-5-ASA from balsalazide';'EX_nchlphncl[fe]','Diet_EX_chlphncl[d]','Nitrosochloramphenicol from chloramphenicol';'EX_bvu[fe]','Diet_EX_srv[d]','(E)-5-(2-Bromovinyl)Uracil from sorivudine';'EX_dopa[fe]','Diet_EX_34dhphe[d]','Dopamine from levodopa';'EX_mtym[fe]','Diet_EX_34dhphe[d]','m-Tyramine from levodopa';'EX_pcresol[fe]','Diet_EX_4hphac[d]','p-cresol from 4-hydroxyphenylacetate';'EX_cholate[fe]','Diet_EX_tchola[d]','Cholic acid from taurocholic acid'};

f=figure;
cols=[];
for i=2:size(fluxes,2)
    subplot(4,4,i-1)
    
    cols(i-1,:)=[rand rand rand];
    data=str2double(fluxes(3:end,i));
    find1=find(strcmp(fluxLabels(:,1),fluxes{1,i}));
    find2=find(strcmp(fluxLabels(:,2),fluxes{2,i}));
    hold on
    violinplot(data, {'Drug conversion potential'},'ViolinColor',cols(i-1,:));
    set(gca, 'FontSize', 10)
    h=title(fluxLabels{intersect(find1,find2),3});
    set(h,'interpreter','none')
    set(gca,'TickLabelInterpreter','none')
end
f.Renderer='painters';

