
%% Compute conversion capability
global CBT_LP_SOLVER
solver = CBT_LP_SOLVER;
% initialize parallel pool
if numWorkers > 0
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end

resultsFolder = [rootDir 'ComputeDrugReactions' filesep];

database=loadVMHDatabase;

WesternDiet = readtable('WesternDietAGORA2.txt', 'Delimiter', '\t');
WesternDiet=table2cell(WesternDiet);
WesternDiet=cellstr(string(WesternDiet));

taxonomy = readtable('AGORA2_infoFile.xlsx', 'ReadVariableNames', false);
taxonomy = table2cell(taxonomy);

for i = 2:size(taxonomy,1)
    model=readCbModel([refinedFolder taxonomy{i,1} '.mat']);
    model.lb(find(strncmp(model.rxns,'sink_',5)))=-1;
    models{i,1}=model;
end

%% compute conversion potential for Figure 3e)

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
        i
        changeCobraSolver(solver, 'LP');
        % prevent creation of log files
        changeCobraSolverParams('LP', 'logFile', 0);
        model=models{i,1};
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
    save('drugPredictions.mat','drugPredictions');
end

% remove the ones not producing anything
cnt=1;
delArray=[];
for j=2:size(drugPredictions,1)
    if abs(sum(cell2mat(drugPredictions(j,2:end))))<0.0001
        delArray(cnt,1)=j;
        cnt=cnt+1;
    end
end
drugPredictions(delArray,:)=[];
cnt=1;
delArray=[];
for j=2:size(drugPredictions,2)
    if abs(sum(cell2mat(drugPredictions(2:end,j))))<0.0001
        delArray(cnt,1)=j;
        cnt=cnt+1;
    end
end
drugPredictions(:,delArray)=[];
cell2csv([resultsFolder 'AGORA2_DrugConversion.csv'],drugPredictions);
% remove glucuronic acid-not needed for plot
drugPredictions(:,2)=[];

taxonomy_reduced=taxonomy;
[C,IA] = setdiff(taxonomy_reduced(:,1),drugPredictions(:,1),'stable');
taxonomy_reduced(IA(2:end),:)=[];
taxonomy_reduced(:,2:4)=strrep(taxonomy_reduced(:,2:4),',','_');
cell2csv([resultsFolder 'AGORA2_Reconstructions_Information_reduced.csv'],taxonomy_reduced);

%% summarize drug conversion capabilities by phylum-Figure 3e

% get organisms with comparative ggitenomics
genomeAnnotation = readtable('gapfilledGenomeAnnotation.txt', 'ReadVariableNames', false, 'Delimiter', 'tab','TreatAsEmpty',['UND. -60001','UND. -2011','UND. -62011']);
genomeAnnotation = table2cell(genomeAnnotation);
orgs=unique(genomeAnnotation(:,1));

drugPredictions = readtable([resultsFolder 'AGORA2_DrugConversion.csv'], 'ReadVariableNames', false);
drugPredictions = table2cell(drugPredictions);
drugPredictions(:,2)=[];
[C,IA] = setdiff(drugPredictions(:,1),orgs,'stable');
drugPredictions(IA(2:end),:)=[];
taxCol=find(strcmp(taxonomy_reduced(1,:),'Phylum'));
get_phyla=unique(taxonomy_reduced(2:end,taxCol));
capabilities={};
for j=2:size(drugPredictions,2)
    capabilities{1,j}=drugPredictions{1,j};
end
capabilities{1,1}='Phylum';
for i=1:length(get_phyla)
    capabilities{i+1,1}=get_phyla{i};
    capabilities(i+1,2:end)={'0'};
    for j=2:size(drugPredictions,2)
        capabilities{i+1,j}=str2double(capabilities{i+1,j});
    end
    alltax=taxonomy_reduced(find(strcmp(taxonomy_reduced(:,taxCol),get_phyla{i})),1);
    for k=1:length(alltax)
        findtax=find(strcmp(drugPredictions(:,1),alltax{k}));
        for j=2:size(drugPredictions,2)
            capabilities{i+1,j}=capabilities{i+1,j} + str2double(drugPredictions{findtax,j});
        end
    end
end

capabilities(1,:)=strrep(capabilities(1,:),'-','_');
capabilities(1,:)=strrep(capabilities(1,:),' ','_');
capabilities=cell2table(capabilities);
writetable(capabilities,[resultsFolder 'AGORA2_Drug_PhylumSummaries'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

%% compute conversion potential for all drugs for comparison against Zimmermann et al. data
resultsFolder = [rootDir filesep 'Drug_metabolism_Comparison_Zimmermann' filesep];

drugExchanges={'InputMet','InputExchange','OutputMet','OutputExchange';'sn38g','EX_sn38g(e)','sn38','EX_sn38(e)';'dfdcytd','EX_dfdcytd(e)','dfduri','EX_dfduri(e)';'fcsn','EX_fcsn(e)','5fura','EX_5fura(e)';'5fura','EX_5fura(e)','dh5fura','EX_dh5fura(e)';'4hphac','EX_4hphac(e)','pcresol','EX_pcresol(e)';'r788','EX_r788(e)','r406','EX_r406(e)';'bzd','EX_bzd(e)','5asa','EX_5asa(e)';'lactl','EX_lactl(e)','gal','EX_gal(e)';'chlphncl','EX_chlphncl(e)','nchlphncl','EX_nchlphncl(e)';'5asa','EX_5asa(e)','ac5asa','EX_ac5asa(e)';'digoxin','EX_digoxin(e)','dihydro_digoxin','EX_dihydro_digoxin(e)';'srv','EX_srv(e)','bvu','EX_bvu(e)';'34dhphe','EX_34dhphe(e)','dopa','EX_dopa(e)';'tchola','EX_tchola(e)','cholate','EX_cholate(e)';'dopa','EX_dopa(e)','mtym','EX_mtym(e)';'1hibupglu_S','EX_1hibupglu_S(e)','1hibup_S','EX_1hibup_S(e)';'1hmdgluc','EX_1hmdgluc(e)','1ohmdz','EX_1ohmdz(e)';'2hatvacidgluc','EX_2hatvacidgluc(e)','2hatvacid','EX_2hatvacid(e)';'2hatvlacgluc','EX_2hatvlacgluc(e)','2hatvlac','EX_2hatvlac(e)';'2hibupglu_S','EX_2hibupglu_S(e)','2hibup_S','EX_2hibup_S(e)';'2oh_cbz_glc','EX_2oh_cbz_glc(e)','2oh_cbz','EX_2oh_cbz(e)';'2oh_mtz_glc','EX_2oh_mtz_glc(e)','2oh_mtz','EX_2oh_mtz(e)';'3hibupglu_S','EX_3hibupglu_S(e)','3hibup_S','EX_3hibup_S(e)';'3oh_cbz_glc','EX_3oh_cbz_glc(e)','3oh_cbz','EX_3oh_cbz(e)';'3oh_dlor_glc','EX_3oh_dlor_glc(e)','3oh_dlor','EX_3oh_dlor(e)';'3oh_mdea_glc','EX_3oh_mdea_glc(e)','3oh_mdea','EX_3oh_mdea(e)';'3oh_mxn_glc','EX_3oh_mxn_glc(e)','3oh_mxn','EX_3oh_mxn(e)';'4dh_tpno_1glc','EX_4dh_tpno_1glc(e)','4dh_tpno','EX_4dh_tpno(e)';'4hmdgluc','EX_4hmdgluc(e)','4ohmdz','EX_4ohmdz(e)';'4oh_dcf_glc','EX_4oh_dcf_glc(e)','4oh_dcf','EX_4oh_dcf(e)';'4oh_kp_glc','EX_4oh_kp_glc(e)','4oh_kp','EX_4oh_kp(e)';'4oh_levole_glc','EX_4oh_levole_glc(e)','4oh_levole','EX_4oh_levole(e)';'4oh_meth_glc','EX_4oh_meth_glc(e)','4oh_meth','EX_4oh_meth(e)';'4oh_propl_glc','EX_4oh_propl_glc(e)','4oh_propl','EX_4oh_propl(e)';'4oh_trz_glc','EX_4oh_trz_glc(e)','4oh_trz','EX_4oh_trz(e)';'4oh_vcz_glc','EX_4oh_vcz_glc(e)','4oh_vcz','EX_4oh_vcz(e)';'dh5fura','EX_dh5fura(e)','5fura','EX_5fura(e)';'5oh_sulfp_glc','EX_5oh_sulfp_glc(e)','5oh_sulfp','EX_5oh_sulfp(e)';'5ohfvsglu','EX_5ohfvsglu(e)','5ohfvs','EX_5ohfvs(e)';'6bhglzglc','EX_6bhglzglc(e)','6bhglz','EX_6bhglz(e)';'6ohfvsglu','EX_6ohfvsglu(e)','6ohfvs','EX_6ohfvs(e)';'7bhglzglc','EX_7bhglzglc(e)','7bhglz','EX_7bhglz(e)';'7oh_efv_glc','EX_7oh_efv_glc(e)','7oh_efv','EX_7oh_efv(e)';'814dioh_efv_glc','EX_814dioh_efv_glc(e)','814dioh_efv','EX_814dioh_efv(e)';'8oh_efv_glc','EX_8oh_efv_glc(e)','8oh_efv','EX_8oh_efv(e)';'ac_amn_b_gly_glc','EX_ac_amn_b_gly_glc(e)','ac_amn_b_gly','EX_ac_amn_b_gly(e)';'acmp_glc','EX_acmp_glc(e)','acmp','EX_acmp(e)';'acmpglu','EX_acmpglu(e)','acmp','EX_acmp(e)';'alpz_4oh_glc','EX_alpz_4oh_glc(e)','alpz_4oh','EX_alpz_4oh(e)';'alpz_aoh_glc','EX_alpz_aoh_glc(e)','alpz_aoh','EX_alpz_aoh(e)';'am14_glc','EX_am14_glc(e)','am14','EX_am14(e)';'am1cglc','EX_am1cglc(e)','am1ccs','EX_am1ccs(e)';'am5_glc','EX_am5_glc(e)','am5','EX_am5(e)';'am6_glc','EX_am6_glc(e)','am6','EX_am6(e)';'amio_c_glc','EX_amio_c_glc(e)','amio_c','EX_amio_c(e)';'amio_glc','EX_amio_glc(e)','amio','EX_amio(e)';'amn_b_gly_glc','EX_amn_b_gly_glc(e)','amn_b_gly','EX_amn_b_gly(e)';'amntd_m6_glc','EX_amntd_m6_glc(e)','amntd_m6','EX_amntd_m6(e)';'atvacylgluc','EX_atvacylgluc(e)','atvlac','EX_atvlac(e)';'atvethgluc','EX_atvethgluc(e)','atvacid','EX_atvacid(e)';'atvlacgluc','EX_atvlacgluc(e)','atvlac','EX_atvlac(e)';'bhpm_glc','EX_bhpm_glc(e)','bhpm','EX_bhpm(e)';'bilr_M10','EX_bilr_M10(e)','bilr_355','EX_bilr_355(e)';'bilr_M12','EX_bilr_M12(e)','bilr_M7','EX_bilr_M7(e)';'bilr_M15','EX_bilr_M15(e)','bilr_M14','EX_bilr_M14(e)';'bilr_M16','EX_bilr_M16(e)','bilr_355','EX_bilr_355(e)';'bsn_glc','EX_bsn_glc(e)','bsn','EX_bsn(e)';'bz_glc','EX_bz_glc(e)','bz','EX_bz(e)';'caribupglu_S','EX_caribupglu_S(e)','caribup_s','EX_caribup_s(e)';'cbz_glc','EX_cbz_glc(e)','cbz','EX_cbz(e)';'cd6168_glc','EX_cd6168_glc(e)','cd6168','EX_cd6168(e)';'chlphncl_glc','EX_chlphncl_glc(e)','chlphncl','EX_chlphncl(e)';'clcxb_c_glc','EX_clcxb_c_glc(e)','clcxb_c','EX_clcxb_c(e)';'clcxb_glc','EX_clcxb_glc(e)','clcxb','EX_clcxb(e)';'clobi_glc','EX_clobi_glc(e)','clobi_c','EX_clobi_c(e)';'cvm1gluc','EX_cvm1gluc(e)','crvsm1','EX_crvsm1(e)';'cvm23gluc','EX_cvm23gluc(e)','crvsm23','EX_crvsm23(e)';'czp','EX_czp(e)','7a_czp','EX_7a_czp(e)';'daa_glc','EX_daa_glc(e)','daa','EX_daa(e)';'dcf_glc','EX_dcf_glc(e)','dcf','EX_dcf(e)';'ddea_glc','EX_ddea_glc(e)','ddea','EX_ddea(e)';'des_astzl_glc','EX_des_astzl_glc(e)','des_astzl','EX_des_astzl(e)';'digoxin_glc','EX_digoxin_glc(e)','digoxin','EX_digoxin(e)';'digitoxin','EX_digitoxin(e)','dihydro_digitoxin','EX_dihydro_digitoxin(e)';'dlb_glc','EX_dlb_glc(e)','dlb','EX_dlb(e)';'dnpz_m11','EX_dnpz_m11(e)','dnpz_6des','EX_dnpz_6des(e)';'dnpz_m12','EX_dnpz_m12(e)','dnpz_5des','EX_dnpz_5des(e)';'dnpz_m13','EX_dnpz_m13(e)','dnpz_m9','EX_dnpz_m9(e)';'dnpz_m14','EX_dnpz_m14(e)','dnpz_m9','EX_dnpz_m9(e)';'doh_etr_glc','EX_doh_etr_glc(e)','doh_etr','EX_doh_etr(e)';'doh_vcz_glc','EX_doh_vcz_glc(e)','doh_vcz','EX_doh_vcz(e)';'dxo_glc','EX_dxo_glc(e)','dxo','EX_dxo(e)';'efv_glc','EX_efv_glc(e)','efv','EX_efv(e)';'eltr_glc','EX_eltr_glc(e)','eltr','EX_eltr(e)';'eltr_m3','EX_eltr_m3(e)','sb_611855','EX_sb_611855(e)';'eltr_m4','EX_eltr_m4(e)','sb_y','EX_sb_y(e)';'eztmb_glc','EX_eztmb_glc(e)','eztmb','EX_eztmb(e)';'fvsgluc','EX_fvsgluc(e)','fvs','EX_fvs(e)';'fvstetglu','EX_fvstetglu(e)','fvstet','EX_fvstet(e)';'glc3meacp','EX_glc3meacp(e)','3meacmp','EX_3meacmp(e)';'gltmn_glc','EX_gltmn_glc(e)','gltmn','EX_gltmn(e)';'gmfl_glc','EX_gmfl_glc(e)','gmfl','EX_gmfl(e)';'gmfl_mI_glc','EX_gmfl_mI_glc(e)','gmfl_mI','EX_gmfl_mI(e)';'gmfl_mII_glc','EX_gmfl_mII_glc(e)','gmfl_mII','EX_gmfl_mII(e)';'gmfl_mIII_glc','EX_gmfl_mIII_glc(e)','gmfl_mIII','EX_gmfl_mIII(e)';'gtacmp','EX_gtacmp(e)','tmacmp','EX_tmacmp(e)';'hst_3_glc','EX_hst_3_glc(e)','hst','EX_hst(e)';'hst_37_diglc','EX_hst_37_diglc(e)','hst_7_glc','EX_hst_7_glc(e)';'hst_3glc_7s','EX_hst_3glc_7s(e)','hst_7_s','EX_hst_7_s(e)';'hst_7_glc','EX_hst_7_glc(e)','hst','EX_hst(e)';'hst_7glc_3s','EX_hst_7glc_3s(e)','hst_3_s','EX_hst_3_s(e)';'ibupgluc','EX_ibupgluc(e)','ibup_S','EX_ibup_S(e)';'imn_glc','EX_imn_glc(e)','imn','EX_imn(e)';'inv_m1','EX_inv_m1(e)','inv','EX_inv(e)';'isnzd','EX_isnzd(e)','acisnzd','EX_acisnzd(e)';'isosorbide_5mn_glc','EX_isosorbide_5mn_glc(e)','isosorbide_5mn','EX_isosorbide_5mn(e)';'kprofen_glc','EX_kprofen_glc(e)','kprofen','EX_kprofen(e)';'lstn1gluc','EX_lstn1gluc(e)','lstn','EX_lstn(e)';'lstnm4','EX_lstnm4(e)','lst4exp','EX_lst4EX_p(e)';'lstnm7','EX_lstnm7(e)','lstn','EX_lstn(e)';'mdz_glc','EX_mdz_glc(e)','mdz','EX_mdz(e)';'mdzglc','EX_mdzglc(e)','mdz','EX_mdz(e)';'miso_glc','EX_miso_glc(e)','miso','EX_miso(e)';'mrphn_3glc','EX_mrphn_3glc(e)','mrphn','EX_mrphn(e)';'mrphn_6glc','EX_mrphn_6glc(e)','mrphn','EX_mrphn(e)';'mtz_glc','EX_mtz_glc(e)','mtz','EX_mtz(e)';'N_oh_phtn_glc','EX_N_oh_phtn_glc(e)','N_oh_phtn','EX_N_oh_phtn(e)';'neopront','EX_neopront(e)','sanilamide','EX_sanilamide(e)';'nsldp_m5_glc','EX_nsldp_m5_glc(e)','nsldp_m5','EX_nsldp_m5(e)';'nverp_glc','EX_nverp_glc(e)','nverp','EX_nverp(e)';'nzp','EX_nzp(e)','anzp','EX_anzp(e)';'odsm_egltmn_glc','EX_odsm_egltmn_glc(e)','odsm_egltmn','EX_odsm_egltmn(e)';'odsm_gltmn_glc','EX_odsm_gltmn_glc(e)','odsm_gltmn','EX_odsm_gltmn(e)';'oh_etr_glc','EX_oh_etr_glc(e)','oh_etr','EX_oh_etr(e)';'oh_pbl_glc','EX_oh_pbl_glc(e)','oh_pbl','EX_oh_pbl(e)';'olsa','EX_olsa(e)','5asa','EX_5asa(e)';'phppa_glc','EX_phppa_glc(e)','phppa','EX_phppa(e)';'phtn_glc','EX_phtn_glc(e)','N_oh_phtn','EX_N_oh_phtn(e)';'prob_glc','EX_prob_glc(e)','prob','EX_prob(e)';'pront','EX_pront(e)','sanilamide','EX_sanilamide(e)';'pront_glc','EX_pront_glc(e)','pront','EX_pront(e)';'propl_glc','EX_propl_glc(e)','propl','EX_propl(e)';'prx_mI_glc','EX_prx_mI_glc(e)','prx_mI','EX_prx_mI(e)';'prx_mII_glc','EX_prx_mII_glc(e)','prx_mII','EX_prx_mII(e)';'ptvstgluc','EX_ptvstgluc(e)','ptvst','EX_ptvst(e)';'pvsgluc','EX_pvsgluc(e)','pvs','EX_pvs(e)';'brv','EX_brv(e)','bvu','EX_bvu(e)';'R_6oh_warf_glc','EX_R_6oh_warf_glc(e)','R_6oh_warf','EX_R_6oh_warf(e)';'R_7oh_warf_glc','EX_R_7oh_warf_glc(e)','R_7oh_warf','EX_R_7oh_warf(e)';'R_8oh_warf_glc','EX_R_8oh_warf_glc(e)','R_8oh_warf','EX_R_8oh_warf(e)';'r406_glc','EX_r406_glc(e)','r406','EX_r406(e)';'r529_glc','EX_r529_glc(e)','r529','EX_r529(e)';'rep_glc','EX_rep_glc(e)','rep','EX_rep(e)';'rpn_104557_cb_glc','EX_rpn_104557_cb_glc(e)','rpn_104557','EX_rpn_104557(e)';'rpn_96990_glc','EX_rpn_96990_glc(e)','rpn_96990','EX_rpn_96990(e)';'rpn_oh_glc','EX_rpn_oh_glc(e)','rpn_oh','EX_rpn_oh(e)';'rsvgluc','EX_rsvgluc(e)','rsv','EX_rsv(e)';'S_4oh_warf_glc','EX_S_4oh_warf_glc(e)','S_4oh_warf','EX_S_4oh_warf(e)';'S_6oh_warf_glc','EX_S_6oh_warf_glc(e)','S_6oh_warf','EX_S_6oh_warf(e)';'sch_488128','EX_sch_488128(e)','eztmb','EX_eztmb(e)';'sch_57871_glc','EX_sch_57871_glc(e)','sch_57871','EX_sch_57871(e)';'sfnd_1689_glc','EX_sfnd_1689_glc(e)','sfnd_1689','EX_sfnd_1689(e)';'sftz_glc','EX_sftz_glc(e)','sftz','EX_sftz(e)';'simvgluc','EX_simvgluc(e)','smvacid','EX_smvacid(e)';'smap_glc','EX_smap_glc(e)','smap','EX_smap(e)';'spz_glc','EX_spz_glc(e)','spz','EX_spz(e)';'spz_sfn_glc','EX_spz_sfn_glc(e)','spz_sfn','EX_spz_sfn(e)';'ssz','EX_ssz(e)','5asa','EX_5asa(e)';'stg_m3','EX_stg_m3(e)','stg','EX_stg(e)';'stg_m4','EX_stg_m4(e)','stg','EX_stg(e)';'tat','EX_tat(e)','tlf_a','EX_tlf_a(e)';'tgz_glc','EX_tgz_glc(e)','tgz','EX_tgz(e)';'tlf_a_m1','EX_tlf_a_m1(e)','tlf_a_2','EX_tlf_a_2(e)';'tlf_a_m2','EX_tlf_a_m2(e)','tlf_a_2a','EX_tlf_a_2a(e)';'tlf_a_m3','EX_tlf_a_m3(e)','tlf_a_3','EX_tlf_a_3(e)';'tlf_a_m4','EX_tlf_a_m4(e)','tlf_a_1a','EX_tlf_a_1a(e)';'tlf_a_m4a','EX_tlf_a_m4a(e)','tlf_a_1a','EX_tlf_a_1a(e)';'tlf_a_m5','EX_tlf_a_m5(e)','tlf_a_1b','EX_tlf_a_1b(e)';'tlf_a_m5a','EX_tlf_a_m5a(e)','tlf_a_1b','EX_tlf_a_1b(e)';'tlf_a_m5b','EX_tlf_a_m5b(e)','tlf_a_m5','EX_tlf_a_m5(e)';'tlf_a_m6','EX_tlf_a_m6(e)','tlf_a_4','EX_tlf_a_4(e)';'tlf_a_m9','EX_tlf_a_m9(e)','tlf_a_1x','EX_tlf_a_1x(e)';'tlms_glc','EX_tlms_glc(e)','tlms','EX_tlms(e)';'tolcp_ac_glc','EX_tolcp_ac_glc(e)','tolcp_ac','EX_tolcp_ac(e)';'tolcp_am_glc','EX_tolcp_am_glc(e)','tolcp_am','EX_tolcp_am(e)';'tolcp_glc','EX_tolcp_glc(e)','tolcp','EX_tolcp(e)';'tpno_1glc_4g','EX_tpno_1glc_4g(e)','tpno_4g','EX_tpno_4g(e)';'tpno_4glc','EX_tpno_4glc(e)','tpnoh','EX_tpnoh(e)';'tsacmgluc','EX_tsacmgluc(e)','thsacmp','EX_thsacmp(e)';'regfnb_glc','EX_regfnb_glc(e)','regfnb','EX_regfnb(e)'};
[C,IA]=setdiff(taxonomy(:,1),orgs);
models(IA,:)=[];
drugPredictions={};
drugPredictions{1,1}='ModelID';
for j = 2:length(drugExchanges)
    drugPredictions{1,j}=drugExchanges{j,1};
    fluxesTmp={};
    % parallelisation for faster computation
    parfor i = 1:length(orgs)
        changeCobraSolver(solver, 'LP');
        % prevent creation of log files
        changeCobraSolverParams('LP', 'logFile', 0);
        model=models{i,1};
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
    for i = 1:length(orgs)
        drugPredictions{i,1}=strrep(orgs{i,1},'.mat','');
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
    save('drugPredictions.mat','drugPredictions');
end

% remove the ones not producing anything
cnt=1;
delArray=[];
for j=2:size(drugPredictions,1)
    if abs(sum(cell2mat(drugPredictions(j,2:end))))<0.0001
        delArray(cnt,1)=j;
        cnt=cnt+1;
    end
end
drugPredictions(delArray,:)=[];
cnt=1;
delArray=[];
for j=2:size(drugPredictions,2)
    if abs(sum(cell2mat(drugPredictions(2:end,j))))<0.0001
        delArray(cnt,1)=j;
        cnt=cnt+1;
    end
end
drugPredictions(:,delArray)=[];
cell2csv([resultsFolder 'AGORA2_DrugConversion.csv'],drugPredictions);
% remove glucuronic acid-not needed for plot
drugPredictions(:,2)=[];

taxonomy_reduced=taxonomy;
[C,IA] = setdiff(taxonomy_reduced(:,1),drugPredictions(:,1),'stable');
taxonomy_reduced(IA(2:end),:)=[];
taxonomy_reduced(:,2:4)=strrep(taxonomy_reduced(:,2:4),',','_');
cell2csv([resultsFolder 'AGORA2_Reconstructions_Information_reduced.csv'],taxonomy_reduced);
