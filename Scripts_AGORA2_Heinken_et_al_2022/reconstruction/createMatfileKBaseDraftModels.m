% adapt original AGORA draft reconstructions to nomenclature to be tested

% AGORA 1.03
inputFolder=[rootDir filesep 'original_AGORA_draft_Models' filesep];
dInfo = dir(inputFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

% rename some modelList to account for changed nomenclature
modelnames={'Old AGORA Model ID','New AGORA Model ID';'Bacillus_thermoamylovorans','Bacillus_thermoamylovorans_1A1';'Bacillus_timonensis_JC401','Bacillus_timonensis_10403023';'Bacteroides_massiliensis_B846dnLKV334','Bacteroides_massiliensis_B84634';'Bifidobacterium_stercoris_ATCC_43183','Bifidobacterium_stercoris_DSM_24849';'Corynebacterium_ureicelerivorans','Corynebacterium_ureicelerivorans_DSM_45051';'Dermacoccus_nishinomiyaensis','Dermacoccus_nishinomiyaensis_DSM_20448';'Lactobacillus_reuteri_CF48_3A_','Lactobacillus_reuteri_CF48_3A';'Lactobacillus_reuteri_MM2_3_','Lactobacillus_reuteri_MM2_3';'Lactobacillus_reuteri_MM4_1A_','Lactobacillus_reuteri_MM4_1A'};

% replace '000000.peg' GPRs with PubSeed IDs for some modelList
replaceGPRs={'MicrobeID','PubSeedID';'Abiotrophia_defectiva_ATCC_49176','592010.4';'Achromobacter_xylosoxidans_A8','762376.5';'Acidaminococcus_fermentans_DSM_20731','591001.3';'Acidaminococcus_intestini_RyC_MR95','568816.4';'Acidaminococcus_sp_D21','563191.3';'Acinetobacter_baumannii_AB0057','480119.5';'Acinetobacter_calcoaceticus_PHEA_2','871585.3';'Acinetobacter_haemolyticus_NIPH_261','1217986.3';'Acinetobacter_johnsonii_SH046','575586.4';'Acinetobacter_junii_SH205','575587.3';'Acinetobacter_lwoffii_WJ10621','1046625.3';'Acinetobacter_pittii_ANC_4052','1217689.3';'Acinetobacter_radioresistens_NIPH_2130','1217674.3';'Actinobacillus_pleuropneumoniae_L20','416269.6';'Actinomyces_cardiffensis_F0333','888050.3';'Actinomyces_graevenitzii_C83','435830.3';'Actinomyces_naeslundii_str_Howell_279','1115803.3';'Actinomyces_odontolyticus_ATCC_17982','411466.7';'Actinomyces_oris_K20','871541.3';'Actinomyces_turicensis_ACS_279_V_Col4','883077.3';'Actinomyces_urogenitalis_DSM_15434','525246.3';'Actinomyces_viscosus_C505','562973.4';'Aerococcus_viridans_ATCC_11563','655812.3';'Aeromonas_caviae_Ae398','878320.4';'Aeromonas_hydrophila_subsp_hydrophila_ATCC_7966','380703.7';'Aeromonas_media_WS','1208104.3';'Aeromonas_veronii_B565','998088.4';'Afipia_birgiae_34632','1197906.5';'Aggregatibacter_aphrophilus_NJ8700','634176.4';'Agrobacterium_tumefaciens_F2','1050720.3';'Akkermansia_muciniphila_ATCC_BAA_835','349741.6';'Alcaligenes_faecalis_subsp_faecalis_NCIB_8687','1156918.4';'Alistipes_finegoldii_DSM_17242','679935.3';'Alistipes_indistinctus_YIT_12060','742725.3';'Alistipes_onderdonkii_DSM_19147','1120974.3';'Alistipes_putredinis_DSM_17216','445970.5';'Alistipes_shahii_WAL_8301','717959.3';'Anaerobaculum_hydrogeniformans_OS1_ATCC_BAA_1850','592015.5';'Anaerococcus_hydrogenalis_DSM_7454','561177.4';'Anaerococcus_prevotii_DSM_20548','525919.6';'Anaerococcus_vaginalis_ATCC_51170','655811.4';'Anaerofustis_stercorihominis_DSM_17244','445971.6';'Anaerostipes_caccae_DSM_14662','411490.6';'Anaerostipes_hadrus_DSM_3319','649757.3';'Anaerostipes_sp_3_2_56FAA','665937.3';'Anaerotruncus_colihominis_DSM_17241','445972.6';'Arcanobacterium_haemolyticum_DSM_20595','644284.3';'Arcobacter_butzleri_RM4018','367737.6';'Atopobium_minutum_10063974','997872.3';'Atopobium_parvulum_DSM_20469','521095.6';'Atopobium_rimae_ATCC_49626','553184.4';'Bacillus_amyloliquefaciens_DSM7','692420.6';'Bacillus_atrophaeus_ATCC_49822_1','743710.3';'Bacillus_cereus_AH187_F4810_72','405534.9';'Bacillus_cereus_G9842','405531.7';'Bacillus_clausii_KSM_K16','66692.6';'Bacillus_endophyticus_2102','1196029.3';'Bacillus_fordii_DSM_16014','1121090.3';'Bacillus_halodurans_C_125','272558.8';'Bacillus_licheniformis_ATCC_14580','279010.13';'Bacillus_megaterium_DSM319','592022.4';'Bacillus_mojavensis_RO_H_1','1051501.3';'Bacillus_mycoides_DSM_2048','526997.3';'Bacillus_nealsonii_AAU1','1202533.3';'Bacillus_pseudofirmus_OF4','398511.4';'Bacillus_pumilus_ATCC_7061','536229.4';'Bacillus_sonorensis_L12','1274524.3';'Bacillus_subtilis_str_168','224308.113';'Bacillus_thuringiensis_serovar_thuringiensis_str_T01001','527025.3';'Bacillus_vallismortis_DV1_F_3','1051502.3';'Bacteroides_caccae_ATCC_43185','411901.7';'Pseudoflavonifractor_capillosus_strain_ATCC_29799','411467.6';'Bacteroides_cellulosilyticus_DSM_14838','537012.5';'Bacteroides_clarus_YIT_12056','762984.3';'Bacteroides_coprocola_M16_DSM_17136','470145.6';'Bacteroides_coprophilus_DSM_18228','547042.5';'Bacteroides_dorei_DSM_17855','483217.6';'Bacteroides_eggerthii_1_2_48FAA','665953.3';'Bacteroides_eggerthii_DSM_20697','483216.6';'Bacteroides_faecis_MAJ27','1077285.3';'Bacteroides_finegoldii_DSM_17565','483215.6';'Bacteroides_fluxus_YIT_12057','763034.3';'Bacteroides_fragilis_3_1_12','457424.5';'Bacteroides_fragilis_638R','817.1';'Bacteroides_fragilis_NCTC_9343','272559.17';'Bacteroides_fragilis_YCH46','295405.11';'Bacteroides_intestinalis_341_DSM_17393','471870.8';'Bacteroides_massiliensis_B84634','1121098.4';'Bacteroides_nordii_CL02T12C05','997884.3';'Bacteroides_oleiciplenus_YIT_12058','742727.4';'Bacteroides_ovatus_ATCC_8483','411476.11';'Bacteroides_ovatus_SD_CC_2a','702444.3';'Bacteroides_ovatus_SD_CMC_3f','702443.3';'Bacteroides_pectinophilus_ATCC_43243','483218.5';'Bacteroides_plebeius_M12_DSM_17135','484018.6';'Bacteroides_salyersiae_WAL_10018','1121101.3';'Bacteroides_sp_1_1_14','469585.3';'Bacteroides_sp_1_1_30','457387.3';'Bacteroides_sp_1_1_6','469586.3';'Bacteroides_sp_2_1_22','469588.3';'Bacteroides_sp_2_1_33B','469589.3';'Bacteroides_sp_2_1_7','457388.5';'Bacteroides_sp_2_2_4','469590.5';'Bacteroides_sp_20_3','469591.4';'Bacteroides_sp_3_1_19','469592.4';'Bacteroides_sp_3_1_23','457390.3';'Bacteroides_sp_3_1_33FAA','457391.3';'Bacteroides_sp_3_1_40A','469593.3';'Bacteroides_sp_3_2_5','457392.3';'Bacteroides_sp_4_1_36','457393.3';'Bacteroides_sp_4_3_47FAA','457394.3';'Bacteroides_sp_9_1_42FAA','457395.6';'Bacteroides_sp_D1','556258.5';'Bacteroides_sp_D2','556259.3';'Bacteroides_sp_D20','585543.3';'Bacteroides_sp_D22','585544.3';'Bacteroides_stercoris_ATCC_43183','449673.7';'Bacteroides_thetaiotaomicron_VPI_5482','226186.12';'Bacteroides_uniformis_ATCC_8492','411479.1';'Bacteroides_ureolyticus_DSM_20703','1121102.3';'Bacteroides_vulgatus_ATCC_8482','435590.9';'Bacteroides_xylanisolvens_SD_CC_1b','702447.3';'Bacteroides_xylanisolvens_XB1A','657309.4';'Barnesiella_intestinihominis_YIT_11860','742726.3';'Bartonella_quintana_Toulouse','283165.4';'Bifidobacterium_adolescentis_ATCC_15703','367928.6';'Bifidobacterium_angulatum_DSM_20098','518635.5';'Bifidobacterium_animalis_lactis_AD011','442563.3';'Bifidobacterium_animalis_lactis_BB_12','552531.3';'Bifidobacterium_animalis_lactis_Bi_07','742729.3';'Bifidobacterium_animalis_lactis_Bl_04_ATCC_SD5219','580050.3';'Bifidobacterium_animalis_lactis_DSM_10140','555970.3';'Bifidobacterium_bifidum_BGN4','484020.3';'Bifidobacterium_bifidum_NCIMB_41171','398513.5';'Bifidobacterium_bifidum_PRL2010','702459.3';'Bifidobacterium_bifidum_S17','883062.5';'Bifidobacterium_breve_HPH0326','1203540.3';'Bifidobacterium_breve_UCC2003_NCIMB8807','326426.4';'Bifidobacterium_catenulatum_DSM_16992','566552.4';'Bifidobacterium_dentium_ATCC_27678','473819.7';'Bifidobacterium_gallicum_DSM_20093','561180.4';'Bifidobacterium_longum_DJO10A','205913.11';'Bifidobacterium_longum_NCC2705','206672.9';'Bifidobacterium_longum_infantis_157F_NC','565040.3';'Bifidobacterium_longum_infantis_ATCC_15697','391904.5';'Bifidobacterium_longum_longum_ATCC_55813','548480.3';'Bifidobacterium_longum_longum_CCUG_52486','537937.5';'Bifidobacterium_longum_longum_BBMN68','890402.3';'Bifidobacterium_longum_longum_JCM_1217','565042.3';'Bifidobacterium_longum_longum_JDM301','759350.3';'Bifidobacterium_pseudocatenulatum_DSM_20438','547043.6';'Bifidobacterium_mongoliense_DSM_21395','547043.8';'Bifidobacterium_thermophilum_RBL67','1254439.12';'Bilophila_wadsworthia_3_1_6','563192.3';'Blautia_hansenii_VPI_C7_24_DSM_20583','537007.6';'Blautia_hydrogenotrophica_DSM_10507','476272.5';'Brachybacterium_paraconglomeratum_LC44','1064537.3';'Brachyspira_pilosicoli_B2904','1133568.3';'Bradyrhizobium_elkanii_USDA_76','398525.3';'Bradyrhizobium_japonicum_USDA_6','1037409.3';'Brevibacillus_agri_BAB_2500','1246477.3';'Brevibacillus_borstelensis_AK1','1300222.3';'Brevibacillus_brevis_NBRC_100599','358681.3';'Brevibacterium_casei_S18','1229781.4';'Brevibacterium_linens_BL2','321955.4';'Brevibacterium_massiliense_5401308','1176165.3';'Brevundimonas_diminuta_470_4','1035191.3';'Marvinbryantia_formatexigens_I_52_DSM_14469','478749.5';'Bulleidia_extructa_W1219','679192.3';'Burkholderia_cepacia_GG4','1009846.3';'Burkholderiales_bacterium_1_1_47','469610.4';'Butyricicoccus_pullicaecorum_1_2','1203606.4';'Butyricimonas_synergistica_DSM_23225','1121129.3';'Butyrivibrio_crossotus_DSM_2876','511680.4';'Butyrivibrio_fibrisolvens_16_4','657324.3';'Campylobacter_coli_JV20','864566.3';'Campylobacter_concisus_13826','360104.9';'Campylobacter_curvus_525_92','360105.8';'Campylobacter_fetus_subsp_fetus_82_40','360106.6';'Campylobacter_gracilis_RM3268','553220.3';'Campylobacter_hominis_ATCC_BAA_381','360107.7';'Campylobacter_jejuni_jejuni_81_176','354242.17';'Campylobacter_jejuni_jejuni_ICDCCJ07004','1295400.3';'Campylobacter_jejuni_jejuni_M1','645464.3';'Campylobacter_jejuni_jejuni_NCTC_11168','192222.6';'Campylobacter_lari_RM2100','306263.5';'Campylobacter_rectus_RM3267','553218.4';'Campylobacter_showae_CSUNSWCD','1244083.3';'Campylobacter_upsaliensis_JV21','888826.3';'Capnocytophaga_granulosa_ATCC_51502','641143.3';'Capnocytophaga_ochracea_DSM_7271','521097.5';'Capnocytophaga_sputigena_ATCC_33612','553177.6';'Catenibacterium_mitsuokai_DSM_15897','451640.5';'Cedecea_davisae_DSM_4568','566551.4';'Citrobacter_freundii_ATCC_8090','1006003.3';'Citrobacter_koseri_ATCC_BAA_895','290338.8';'Citrobacter_sp_30_2','469595.3';'Citrobacter_youngae_ATCC_29220','500640.5';'Clostridiales_sp_1_7_47FAA','457421.5';'Clostridium_acetobutylicum_ATCC_824','272562.8';'Clostridium_asparagiforme_DSM_15981','518636.5';'Clostridium_bartlettii_DSM_16795','445973.7';'Clostridium_beijerinckii_NCIMB_8052','290402.41';'Clostridium_bolteae_ATCC_BAA_613','411902.9';'Clostridium_botulinum_A_str_ATCC_19397','441770.6';'Clostridium_botulinum_A_str_ATCC_3502','413999.7';'Clostridium_botulinum_A_str_Hall','441771.6';'Clostridium_botulinum_Bf','445336.4';'Clostridium_botulinum_E1_str_BoNT_E_Beluga','536233.3';'Clostridium_botulinum_F_str_230613','758678.3';'Clostridium_butyricum_DSM_10702','1316931.3';'Clostridium_butyricum_E4_str_BoNT_E_BL5262','632245.3';'Clostridium_celatum_DSM_1785','545697.3';'Clostridium_citroniae_WAL_17108','742733.3';'Clostridium_clariflavum_DSM_19732','720554.3';'Clostridium_clostridioforme_CM201','999410.3';'Clostridioides_difficile_CD196','645462.3';'Clostridioides_difficile_NAP07','525258.3';'Clostridioides_difficile_NAP08','525259.3';'Clostridioides_difficile_R20291','645463.3';'Clostridium_glycolicum_ATCC_14880','1121315.3';'Clostridium_hathewayi_12489931','999412.4';'Clostridium_hiranonis_TO_931_DSM_13275','500633.7';'Clostridium_hylemonae_DSM_15053','553973.6';'Clostridium_innocuum_2959','999413.4';'Clostridium_leptum_DSM_753','428125.8';'Clostridium_methylpentosum_R2_DSM_5476','537013.3';'Clostridium_nexile_DSM_1787','500632.7';'Clostridium_perfringens_ATCC_13124','195103.1';'Clostridium_ramosum_VPI_0427_DSM_1402','445974.6';'Clostridium_saccharoperbutylacetonicum_N1_4_HMT','931276.5';'Clostridium_sartagoforme_AAU1','1202534.3';'Clostridium_scindens_ATCC_35704','411468.9';'Clostridium_sp_7_2_43FAA','457396.3';'Clostridium_sp_L2_50','411489.7';'Clostridium_sp_M62_1','411486.3';'Clostridium_sp_SS2_1','411484.7';'Clostridium_sp_SY8519','1042156.4';'Clostridium_spiroforme_DSM_1552','428126.7';'Clostridium_sporogenes_ATCC_15579','471871.7';'Clostridium_sporosphaeroides_DSM_1294','1121334.3';'Clostridium_sticklandii_DSM_519','499177.3';'Clostridium_symbiosum_WAL_14163','742740.3';'Clostridium_symbiosum_WAL_14673','742741.3';'Clostridium_tyrobutyricum_DSM_2637','1121342.4';'Collinsella_aerofaciens_ATCC_25986','411903.6';'Collinsella_intestinalis_DSM_13280','521003.7';'Collinsella_stercoris_DSM_13279','445975.6';'Collinsella_tanakaei_YIT_12063','742742.3';'Comamonas_testosteroni_CNB_2','688245.4';'Coprobacillus_cateniformis_29_1','469596.3';'Coprococcus_catus_GD_7','717962.3';'Coprococcus_comes_ATCC_27758','470146.3';'Coprococcus_eutactus_ATCC_27759','411474.6';'Corynebacterium_ammoniagenes_DSM_20306','649754.3';'Corynebacterium_amycolatum_SK46','553204.6';'Corynebacterium_aurimucosum_ATCC_700975','548476.7';'Corynebacterium_durum_F0235','1035195.3';'Corynebacterium_glucuronolyticum_ATCC_51867','548477.4';'Corynebacterium_kroppenstedtii_DSM_44385','645127.4';'Corynebacterium_propinquum_DSM_44285','1121367.3';'Corynebacterium_striatum_ATCC_6940','525268.3';'Corynebacterium_tuberculostearicum_SK141','553206.4';'Corynebacterium_ulcerans_809','945711.3';'Cronobacter_sakazakii_ATCC_BAA_894','290339.8';'Cryptobacterium_curtum_DSM_15641','469378.5';'Curtobacterium_flaccumfaciens_UCD_AKU','1292022.3';'Delftia_acidovorans_SPH_1','398578.5';'Desulfovibrio_piger_ATCC_29098','411464.8';'Desulfovibrio_sp_3_1_syn3','457398.5';'Dialister_invisus_DSM_15470','592028.3';'Dialister_succinatiphilus_YIT_11850','742743.3';'Dietzia_cinnamea_P4','910954.3';'Dorea_formicigenerans_ATCC_27755','411461.4';'Dorea_longicatena_DSM_13814','411462.6';'Dyadobacter_beijingensis_DSM_21582','1121482.3';'Dyadobacter_fermentans_DSM_18053','471854.5';'Dysgonomonas_gadei_ATCC_BAA_286','742766.3';'Edwardsiella_tarda_ATCC_23685','500638.3';'Edwardsiella_tarda_FL6_60','718251.5';'Eggerthella_lenta_DSM_2243','479437.5';'Eggerthella_sp_1_3_56FAA','665943.3';'Eggerthella_sp_YY7918','502558.4';'Eggerthia_catenaformis_OT_569','999415.3';'Eikenella_corrodens_ATCC_23834','546274.4';'Enterobacter_aerogenes_KCTC_2190','1028307.3';'Enterobacter_asburiae_LF7a','640513.3';'Enterobacter_cancerogenus_ATCC_35316','500639.8';'Enterobacter_cloacae_EcWSU1','1045856.3';'Enterobacter_hormaechei_YT2','1259823.3';'Enterobacter_hormaechei_YT3','1260282.3';'Enterobacteriaceae_bacterium_9_2_54FAA','469613.3';'Enterococcus_asini_ATCC_700915','1158606.4';'Enterococcus_avium_ATCC_14025','1140002.4';'Enterococcus_caccae_ATCC_BAA_1240','1158612.4';'Enterococcus_casseliflavus_ATCC_12755','888066.3';'Enterococcus_cecorum_DSM_20682','1121864.4';'Enterococcus_dispar_ATCC_51266','1139219.4';'Enterococcus_durans_ATCC_6056','1140001.4';'Enterococcus_faecalis_OG1RF_ATCC_47077','474186.5';'Enterococcus_faecalis_TX0104','491074.3';'Enterococcus_faecalis_TX1322','525278.3';'Enterococcus_faecalis_TX2134','749518.3';'Enterococcus_faecalis_V583','226185.9';'Enterococcus_faecium_TX1330','525279.3';'Enterococcus_gallinarum_EG2','565653.4';'Enterococcus_hirae_ATCC_9790','768486.3';'Enterococcus_phoeniculicola_ATCC_BAA_412','1158610.4';'Enterococcus_sp_7L76','657310.3';'Erysipelotrichaceae_bacterium_sp_3_1_53','658659.3';'Erysipelotrichaceae_bacterium_5_2_54FAA','552396.3';'Escherichia_albertii_TW07627','502347.3';'Escherichia_coli_O157_H7_str_Sakai','386585.9';'Escherichia_coli_SE11','409438.11';'Escherichia_coli_str_K_12_substr_MG1655','511145.12';'Escherichia_coli_UTI89_UPEC','364106.8';'Escherichia_fergusonii_ATCC_35469','585054.5';'Escherichia_hermannii_NBRC_105704','1115512.3';'Escherichia_sp_1_1_43','457400.3';'Escherichia_sp_3_2_53FAA','469598.5';'Escherichia_sp_4_1_40B','457401.3';'Eubacterium_biforme_DSM_3989','518637.5';'Eubacterium_cellulosolvens_6','633697.3';'Eubacterium_cylindroides_T2_87','717960.3';'Eubacterium_dolichum_DSM_3991','428127.7';'Eubacterium_eligens_ATCC_27750','515620.4';'Eubacterium_hallii_DSM_3353','411469.3';'Eubacterium_callanderi_KIST612','903814.3';'Eubacterium_rectale_ATCC_33656','515619.6';'Eubacterium_rectale_M104_1','657317.3';'Eubacterium_saphenum_ATCC_49989','592031.3';'Eubacterium_siraeum_70_3','657319.3';'Eubacterium_siraeum_DSM_15702','428128.7';'Eubacterium_ventriosum_ATCC_27560','411463.4';'Faecalibacterium_cf_prausnitzii_KLE1255','748224.3';'Faecalibacterium_prausnitzii_A2_165','411483.3';'Faecalibacterium_prausnitzii_L2_6','718252.3';'Faecalibacterium_prausnitzii_M21_2','411485.1';'Faecalibacterium_prausnitzii_SL3_3','657322.3';'Filifactor_alocis_ATCC_35896','546269.5';'Finegoldia_magna_ATCC_29328','334413.6';'Flavonifractor_plautii_ATCC_29863','411475.3';'Fusobacterium_gonidiaformans_3_1_5R','469605.3';'Fusobacterium_gonidiaformans_ATCC_25563','469615.3';'Fusobacterium_mortiferum_ATCC_9817','469616.3';'Fusobacterium_necrophorum_D12','556263.3';'Fusobacterium_nucleatum_subsp_animalis_4_8','469607.3';'Fusobacterium_nucleatum_subsp_animalis_7_1','457405.3';'Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586','190304.8';'Fusobacterium_nucleatum_subsp_vincentii_3_1_36A2','469604.7';'Fusobacterium_nucleatum_subsp_vincentii_4_1_13','469606.3';'Fusobacterium_periodonticum_1_1_41FAA','469621.3';'Fusobacterium_periodonticum_2_1_31','469599.8';'Fusobacterium_russii_ATCC_25533','1278306.3';'Fusobacterium_nucleatum_subsp_animalis_3_1_33','469603.3';'Fusobacterium_nucleatum_subsp_vincentii_3_1_27','469602.3';'Fusobacterium_nucleatum_subsp_animalis_D11','556264.3';'Fusobacterium_ulcerans_ATCC_49185','469617.3';'Fusobacterium_varium_ATCC_27725','469618.3';'Gemella_haemolysans_ATCC_10379','546270.5';'Gemella_morbillorum_M424','562982.3';'Gemella_sanguinis_M325','562983.3';'Gordonia_rubripertincta_NBRC_101908','1077975.4';'Gordonia_terrae_C_6','1316928.3';'Gordonia_terrae_NBRC_100016','1089454.3';'Gordonibacter_pamelaeae_7_10_1_bT_DSM_19378','657308.3';'Granulicatella_adiacens_ATCC_49175','638301.3';'Granulicatella_elegans_ATCC_700633','626369.3';'Grimontia_hollisae_CIP_101886','675812.3';'Haemophilus_haemolyticus_M19501','1028803.3';'Haemophilus_influenzae_R2846','262727.7';'Haemophilus_parainfluenzae_T3T1','862965.3';'Haemophilus_sputorum_CCUG_13788','1035839.4';'Hafnia_alvei_ATCC_51873','1002364.3';'Helicobacter_bilis_ATCC_43879','613026.4';'Helicobacter_canadensis_MIT_98_5491','537970.9';'Helicobacter_cinaedi_CCUG_18818','537971.5';'Helicobacter_pullorum_MIT_98_5489','537972.5';'Helicobacter_pylori_26695','85962.8';'Helicobacter_winghamensis_ATCC_BAA_430','556267.4';'Holdemania_filiformis_VPI_J1_31B_1_DSM_12042','545696.5';'Kingella_oralis_ATCC_51147','629741.3';'Klebsiella_oxytoca_KCTC_1686','1006551.4';'Klebsiella_pneumoniae_pneumoniae_MGH78578','272620.3';'Klebsiella_variicola_1_1_55','469608.3';'Kocuria_palustris_PEL','1236550.3';'Kocuria_rhizophila_DC2201','378753.5';'Kurthia_massiliensis_JC30','1033739.3';'Kytococcus_sedentarius_DSM_20547','478801.5';'Lachnospiraceae_bacterium_sp_5_1_63FAA','658089.3';'Lachnospiraceae_bacterium_sp_8_1_57FAA','665951.3';'Lactobacillus_acidophilus_ATCC_4796','525306.3';'Lactobacillus_acidophilus_NCFM','272621.13';'Lactobacillus_amylolyticus_DSM_11664','585524.3';'Lactobacillus_amylovorus_GRL_1112','695560.3';'Lactobacillus_animalis_KCTC_3501','930942.3';'Lactobacillus_antri_DSM_16041','525309.3';'Lactobacillus_brevis_ATCC_367','387344.15';'Lactobacillus_brevis_subsp_gravesensis_ATCC_27305','525310.3';'Lactobacillus_buchneri_ATCC_11577','525318.3';'Lactobacillus_casei_ATCC_334','321967.8';'Lactobacillus_casei_casei_BL23','543734.4';'Lactobacillus_coleohominis_101_4_CHN','575594.3';'Lactobacillus_coryniformis_subsp_coryniformis_CECT_5711','1185325.3';'Lactobacillus_crispatus_125_2_CHN','575595.3';'Lactobacillus_curvatus_CRL_705','1074451.3';'Lactobacillus_delbrueckii_subsp_bulgaricus_ATCC_11842','390333.7';'Lactobacillus_delbrueckii_subsp_bulgaricus_ATCC_BAA_365','321956.5';'Lactobacillus_fermentum_ATCC_14931','525325.3';'Lactobacillus_fermentum_IFO_3956','334390.5';'Lactobacillus_gasseri_ATCC_33323','324831.13';'Lactobacillus_gastricus_PS3','1144300.3';'Lactobacillus_helveticus_DPC_4571','405566.6';'Lactobacillus_helveticus_DSM_20075','585520.4';'Lactobacillus_hilgardii_ATCC_8290','525327.3';'Lactobacillus_iners_DSM_13335','525328.4';'Lactobacillus_jensenii_1153','440497.1';'Lactobacillus_johnsonii_DPC_6026','909954.3';'Lactobacillus_johnsonii_NCC_533','257314.6';'Lactobacillus_mucosae_LM1','1130798.3';'Lactobacillus_oris_F0423','944562.4';'Lactobacillus_paracasei_subsp_paracasei_8700_2','537973.8';'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302','525337.3';'Lactobacillus_pentosus_KCA1','1136177.4';'Lactobacillus_plantarum_JDM1','644042.3';'Lactobacillus_plantarum_subsp_plantarum_ATCC_14917','525338.3';'Lactobacillus_plantarum_WCFS1','220668.9';'Lactobacillus_reuteri_CF48_3A','525341.3';'Lactobacillus_reuteri_F275_JCM_1112','557433.3';'Lactobacillus_reuteri_MM2_3','585517.3';'Lactobacillus_reuteri_MM4_1A','548485.3';'Lactobacillus_reuteri_SD2112_ATCC_55730','491077.3';'Lactobacillus_rhamnosus_GG_ATCC_53103','568703.3';'Lactobacillus_rhamnosus_LMS2_1','525361.3';'Lactobacillus_ruminis_ATCC_25644','525362.3';'Lactobacillus_sakei_subsp_sakei_23K','314315.12';'Lactobacillus_salivarius_salivarius_UCC118','362948.14';'Lactobacillus_ultunensis_DSM_16047','525365.3';'Lactobacillus_vaginalis_ATCC_49540','525366.3';'Lactococcus_garvieae_ATCC_49156','420889.6';'Lactococcus_lactis_subsp_lactis_Il1403','272623.7';'Lactococcus_raffinolactis_4877','1215915.3';'Laribacter_hongkongensis_HLHK9','557598.3';'Lautropia_mirabilis_ATCC_51599','887898.3';'Leptotrichia_buccalis_C_1013_b','523794.5';'Leuconostoc_argentinum_KCTC_3773','886872.3';'Leuconostoc_gelidum_JB7','1229756.3';'Leuconostoc_mesenteroides_subsp_cremoris_ATCC_19254','586220.3';'Listeria_grayi_DSM_20601','525367.9';'Listeria_monocytogenes_Finland_1988','393127.4';'Listeria_monocytogenes_FSL_R2_561','393126.4';'Listeria_monocytogenes_J0161_FSL_R2_499','393130.8';'Listeria_monocytogenes_serotype_7_SLCC_2482','1639.166';'Listeria_monocytogenes_SLCC_2378','879088.3';'Listeria_monocytogenes_SLCC_2540','879089.3';'Listeria_monocytogenes_SLCC_7179','879090.3';'Listeria_monocytogenes_4b_F2365','265669.9';'Lysinibacillus_fusiformis_ZB2','1231627.3';'Lysinibacillus_sphaericus_C3_41','444177.5';'Megamonas_funiformis_YIT_11815','742816.3';'Megamonas_hypermegale_ART12_1','657316.3';'Megasphaera_elsdenii_DSM_20460','907.4';'Mesorhizobium_loti_MAFF303099','266835.9';'Methanobrevibacter_smithii_ATCC_35061','420247.6';'Methanosphaera_stadtmanae_DSM_3091','339860.6';'Methylobacterium_mesophilicum_SR1_6_6','908290.3';'Methylobacterium_populi_BJ001','441620.7';'Methylobacterium_radiotolerans_JCM_2831','426355.14';'Methyloversatilis_universalis_FAM5','1000565.3';'Microbacterium_paraoxydans_77MFTsu3_2','1151126.3';'Micrococcus_luteus_NCTC_2665','465515.4';'Micromonospora_aurantiaca_ATCC_27029','644283.3';'Mitsuokella_multacida_DSM_20544','500635.8';'Mobiluncus_curtisii_ATCC_43063','548479.6';'Coprobacillus_sp_D7','556270.3';'Moraxella_catarrhalis_RH4','749219.3';'Morganella_morganii_subsp_morganii_KT','1124991.3';'Mycobacterium_avium_subsp_avium_ATCC_25291','553481.3';'Mycobacterium_fortuitum_subsp_fortuitum_DSM_46621','1214102.3';'Mycoplasma_hominis_ATCC_23114','347256.5';'Mycoplasma_pneumoniae_309','1112856.4';'Neisseria_cinerea_ATCC_14685','546262.4';'Neisseria_elongata_subsp_glycolytica_ATCC_29315','546263.3';'Neisseria_flavescens_SK114','596320.3';'Neisseria_macacae_ATCC_33926','997348.4';'Neisseria_mucosa_ATCC_25996','546266.6';'Neisseria_subflava_NJ9703','546268.4';'Ochrobactrum_anthropi_ATCC_49188','439375.7';'Ochrobactrum_intermedium_LMG_3301','641118.3';'Odoribacter_laneus_YIT_12061','742817.3';'Odoribacter_splanchnicus_1651_6_DSM_20712','709991.3';'Olsenella_uli_DSM_7084','633147.4';'Oribacterium_sinus_F0268','585501.3';'Oxalobacter_formigenes_HOxBLS','556268.6';'Oxalobacter_formigenes_OXCC13','556269.4';'Paenibacillus_alvei_DSM_29','1206781.3';'Paenibacillus_barengoltzii_G22','1235795.3';'Paenibacillus_daejeonensis_DSM_15491','1122917.3';'Paenibacillus_lactis_154','743719.3';'Pantoea_agglomerans_IG1','1110694.4';'Parabacteroides_distasonis_ATCC_8503','435591.13';'Parabacteroides_goldsteinii_dnLKV18','1235789.3';'Parabacteroides_johnsonii_DSM_18315','537006.5';'Parabacteroides_merdae_ATCC_43184','411477.4';'Parabacteroides_sp_D13','563193.3';'Paraprevotella_clara_YIT_11840','762968.3';'Paraprevotella_xylaniphila_YIT_11841','762982.3';'Parasutterella_excrementihominis_YIT_11859','762966.3';'Parvimonas_micra_ATCC_33270','411465.1';'Pediococcus_acidilactici_7_4','563194.3';'Pediococcus_acidilactici_DSM_20284','862514.3';'Pediococcus_pentosaceus_ATCC_25745','278197.12';'Peptoniphilus_harei_ACS_146_V_Sch2b','908338.3';'Peptoniphilus_indolicus_ATCC_29427','997350.3';'Peptoniphilus_lacrimalis_DSM_7455','1122949.3';'Peptoniphilus_timonensis_JC401','1095770.3';'Peptostreptococcus_anaerobius_DSM_2949','1122950.3';'Peptostreptococcus_stomatis_DSM_17678','596315.3';'Phascolarctobacterium_succinatutens_YIT_12067','626939.3';'Plesiomonas_shigelloides_302_73','1315976.3';'Porphyromonas_asaccharolytica_DSM_20707','879243.3';'Porphyromonas_endodontalis_ATCC_35406','553175.3';'Porphyromonas_gingivalis_W83','242619.8';'Porphyromonas_somerae_DSM_23386','1122975.3';'Porphyromonas_uenonis_60_3','596327.3';'Prevotella_bivia_DSM_20514','868129.3';'Prevotella_bryantii_B14','752555.5';'Prevotella_buccae_ATCC_33574','873513.3';'Prevotella_conceptionensis_9403948','1197728.3';'Prevotella_copri_CB7_DSM_18205','537011.5';'Prevotella_denticola_F0289','767031.3';'Prevotella_disiens_FB035_09AN','866771.4';'Prevotella_intermedia_17','246198.1';'Prevotella_loescheii_DSM_19665','1122985.3';'Prevotella_melaninogenica_ATCC_25845','553174.6';'Prevotella_nanceiensis_DSM_19126','1122988.3';'Prevotella_nigrescens_ATCC_33563','997352.4';'Prevotella_oralis_ATCC_33269','873533.3';'Prevotella_pallens_ATCC_700821','997353.4';'Prevotella_ruminicola_23','264731.4';'Prevotella_salivae_DSM_15606','888832.3';'Prevotella_stercorea_DSM_18206','1002367.3';'Prevotella_timonensis_CRIS_5C_B1','679189.3';'Prevotella_veroralis_DSM_19559','1122993.3';'Propionibacterium_acidipropionici_ATCC_4875','1171373.8';'Cutibacterium_acnes_KPA171202','267747.3';'Propionibacterium_avidum_44067','1170318.3';'Propionibacterium_freudenreichii_subsp_shermanii_CIRM_BIA1','754252.3';'Propionibacterium_propionicum_F0230a','767029.3';'Proteus_mirabilis_ATCC_29906','525369.4';'Proteus_penneri_ATCC_35198','471881.3';'Providencia_alcalifaciens_DSM_30120','520999.6';'Providencia_rettgeri_DSM_1131','521000.6';'Providencia_rustigianii_DSM_4541','500637.6';'Providencia_stuartii_ATCC_25827','471874.6';'Pseudomonas_aeruginosa_NCGM2_S1','1089456.3';'Pseudomonas_alcaliphila_34','1248438.3';'Pseudomonas_fluorescens_PfO_1','205922.5';'Pseudomonas_monteilii_QM','1123524.3';'Pseudomonas_putida_F1','351746.6';'Pseudomonas_stutzeri_DSM_4166','996285.3';'Pseudoramibacter_alactolyticus_ATCC_23263','887929.3';'Pyramidobacter_piscolens_W5455','352165.3';'Ralstonia_pickettii_5_7_47FAA','658664.3';'Rhodococcus_equi_ATCC_33707','525370.3';'Rhodococcus_erythropolis_PR4','234621.6';'Roseburia_hominis_A2_183','585394.18';'Roseburia_intestinalis_L1_82','536231.5';'Roseburia_inulinivorans_DSM_16841','622312.4';'Rothia_aeria_F0474','1125724.3';'Rothia_dentocariosa_ATCC_17931','762948.4';'Rothia_mucilaginosa_DY_18','680646.3';'Rudanella_lutea_DSM_19387','1089547.3';'Ruminococcaceae_bacterium_D16','552398.3';'Ruminococcus_albus_7','697329.1';'Ruminococcus_bromii_L2_63','657321.5';'Ruminococcus_flavefaciens_FD_1','641112.4';'Blautia_gnavus_ATCC_29149','411470.6';'Ruminococcus_lactaris_ATCC_29176','471875.6';'Blautia_obeum_A2_162','657314.3';'Blautia_obeum_ATCC_29174','411459.7';'Ruminococcus_champanellensis_18P13','213810.4';'Ruminococcus_sp_5_1_39BFAA','457412.4';'Ruminococcus_sp_SR1_5','657323.3';'Blautia_torques_ATCC_27756','411460.6';'Blautia_torques_L2_14','657313.3';'Salmonella_enterica_enterica_sv_Typhimurium_LT2','99287.12';'Scardovia_inopinata_F0304','641146.3';'Schlesneria_paludicola_DSM_18645','1123242.3';'Serratia_marcescens_subsp_marcescens_Db11','615.1';'Shigella_dysenteriae_Sd197','300267.13';'Shigella_flexneri_2002017','591020.3';'Shigella_sonnei_Ss046','300269.12';'Slackia_exigua_ATCC_700122','649764.3';'Slackia_piriformis_YIT_12062','742818.3';'Solobacterium_moorei_F0204','706433.3';'Spirosoma_linguale_DSM_74','504472.7';'Staphylococcus_arlettae_CVD059','1212545.3';'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757','451515.3';'Staphylococcus_capitis_QN1','1189311.3';'Staphylococcus_caprae_C87','435838.3';'Staphylococcus_epidermidis_ATCC_12228','176280.1';'Staphylococcus_equorum_subsp_equorum_Mu2','1159488.5';'Staphylococcus_haemolyticus_JCSC1435','279808.8';'Staphylococcus_hominis_subsp_hominis_C80','435837.3';'Staphylococcus_lugdunensis_HKU09_01','698737.3';'Staphylococcus_pettenkoferi_VCU012','904314.5';'Staphylococcus_saprophyticus_subsp_saprophyticus_ATCC_15305','342451.11';'Staphylococcus_simulans_ACS_120_V_Sch1','883166.5';'Staphylococcus_vitulinus_F1028','1167632.5';'Staphylococcus_warneri_SG1','1194526.3';'Staphylococcus_xylosus_NJ','1262650.3';'Stenotrophomonas_maltophilia_D457','1163399.5';'Streptococcus_agalactiae_A909','205921.4';'Streptococcus_anginosus_1_2_62CV','742820.3';'Streptococcus_australis_ATCC_700641','888833.3';'Streptococcus_constellatus_subsp_constellatus_SK53','1095730.3';'Streptococcus_constellatus_subsp_pharyngis_SK1060','1035184.4';'Streptococcus_cristatus_ATCC_51100','889201.3';'Streptococcus_dysgalactiae_subsp_dysgalactiae_ATCC_27957','663952.3';'Streptococcus_equi_subsp_equi_4047','553482.3';'Streptococcus_equinus_ATCC_9812','525379.3';'Streptococcus_gallolyticus_subsp_gallolyticus_ATCC_43143','981539.3';'Streptococcus_gordonii_str_Challis_substr_CH1','467705.9';'Streptococcus_infantarius_subsp_infantarius_ATCC_BAA_102','471872.6';'Streptococcus_infantis_ATCC_700779','889204.3';'Streptococcus_intermedius_JTH08','591365.3';'Streptococcus_mitis_NCTC_12261','246201.6';'Streptococcus_mutans_ATCC_25175','1257041.3';'Streptococcus_oralis_Uo5','927666.3';'Streptococcus_parasanguinis_ATCC_903','888048.3';'Streptococcus_parauberis_KCTC_11537','936154.3';'Streptococcus_peroris_ATCC_700780','888746.3';'Streptococcus_pneumoniae_G54','512566.5';'Streptococcus_pseudopneumoniae_IS7493','1054460.4';'Streptococcus_pyogenes_MGAS9429','370551.4';'Streptococcus_salivarius_JIM8777','347253.5';'Streptococcus_sanguinis_SK36','388919.9';'Streptococcus_sp_2_1_36FAA','469609.3';'Streptococcus_thermophilus_LMG_18311','264199.4';'Streptococcus_thoraltensis_DSM_12221','1123318.3';'Streptococcus_uberis_0140J','218495.5';'Streptococcus_vestibularis_F0396','904306.3';'Subdoligranulum_variabile_DSM_15176','411471.5';'Succinatimonas_hippei_YIT_12066','762983.3';'Saccharolobus_solfataricus_P2','273057.1';'Sutterella_parvirubra_YIT_11816','762967.3';'Sutterella_wadsworthensis_3_1_45B','742821.3';'Tannerella_forsythia_ATCC_43037','203275.8';'Tropheryma_whipplei_str_Twist','203267.6';'Ureaplasma_parvum_serovar_1_str_ATCC_27813','515608.4';'Ureaplasma_urealyticum_serovar_8_str_ATCC_27618','626095.4';'Ureibacillus_thermosphaericus_str_Thermo_BF','1160707.3';'Veillonella_atypica_ACS_049_V_Sch6','866776.4';'Veillonella_dispar_ATCC_17748','546273.3';'Veillonella_parvula_Te3_DSM_2008','479436.6';'Veillonella_ratti_ACS_216_V_Col6b','883156.3';'Veillonella_sp_3_1_44','457416.3';'Veillonella_sp_6_1_27','450749.3';'Vibrio_furnissii_NCTC_11218','903510.3';'Vibrio_mimicus_MB_451','675806.3';'Vibrio_parahaemolyticus_RIMD_2210633','223926.6';'Weissella_cibaria_KACC_11862','911104.3';'Weissella_confusa_LBAE_C39_2','1127131.3';'Weissella_paramesenteroides_ATCC_33313','585506.3';'Yersinia_bercovieri_ATCC_43970','349968.5';'Yersinia_enterocolitica_subsp_enterocolitica_8081','393305.7';'Yersinia_frederiksenii_ATCC_33641','349966.5';'Yersinia_kristensenii_ATCC_33638','527012.3';'Yersinia_pseudotuberculosis_YPIII','502800.6';'Yersinia_rohdei_ATCC_43380','527004.3';'Yokenella_regensburgei_ATCC_43003','1002368.3'};

for i=1:length(modelList)
    % two modelList removed and one replaced
    if ~any(strcmp(modelList{i},{'Staphylococcus_intermedius_ATCC_27335.mat','Paenibacillus_graminis_C4D1M'}))
        try
            model=readCbModel([inputFolder modelList{i,1}]);
        catch
            load([inputFolder modelList{i,1}]);
        end
        outputName=modelList{i,1};
        for j=2:length(modelnames)
            outputName=strrep(outputName,modelnames{j,1},modelnames{j,2});
        end
        findModelID=find(strcmp(replaceGPRs(:,1),strrep(outputName,'.mat','')));
        if ~isempty(findModelID)

            for j=1:length(model.genes)
                model.genes{j}=strrep(model.genes{j},'0000000.0.peg',[replaceGPRs{findModelID,2} '.peg']);
            end
            for j=1:length(model.grRules)
                model.grRules{j}=strrep(model.grRules{j},'0000000.0.peg',[replaceGPRs{findModelID,2} '.peg']);
            end
        end
        
        for j=1:length(model.grRules)
            model.grRules{j}=strrep(model.grRules{j},' and  and  and  and  and  and  and  and','');
            model.grRules{j}=strrep(model.grRules{j},' and )',')');
            model.grRules{j}=strrep(model.grRules{j},' and)',')');
            model.grRules{j}=strrep(model.grRules{j},'and and','and');
            model.grRules{j}=strrep(model.grRules{j},'and  and','and'); 
            model.grRules{j}=strrep(model.grRules{j},'andandandand )',')'); 
            model.grRules{j}=strrep(model.grRules{j},'andandandand','and'); 
            model.grRules{j}=strrep(model.grRules{j},'andand','and');
            model.grRules{j}=strrep(model.grRules{j},'andand','and');
            model.grRules{j}=strrep(model.grRules{j},' and )',')');
        end
        
        % repair gene rules
        if isfield(model,'rules')
            model=rmfield(model,'rules');
        end
        
        % rebuild rules field
        model = generateRules(model);
        
        model = translateDraftReconstruction(model);
        
        errs=verifyModel(model);
        finderrs=fieldnames(errs);
        if ~isempty(finderrs)
            error('problem')
        end
        writeCbModel(model,'fileName',[translatedDraftdsFolder outputName],'format','mat');
    end
end
