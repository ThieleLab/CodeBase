function [modelOut]=useWesternDiet_AGORA(modelIn)
% assign a Western diet for the AGORA microbes
% used as constraints to predict the growth rates shown in Table S6 in
% Magnusdottir et al., Nat Biotech 2017
% defined similarly as in Heinken and Thiele, AEM 2015
% anaerobic medium

model=modelIn;
model=changeRxnBounds(model,model.rxns(strmatch('EX_',model.rxns)),0,'l');
% simple sugars and starch
model=changeRxnBounds(model,'EX_fru(e)',-0.148986,'l');
model=changeRxnBounds(model,'EX_glc(e)',-0.148986,'l');
model=changeRxnBounds(model,'EX_gal(e)',-0.148986,'l');
model=changeRxnBounds(model,'EX_man(e)',-0.148986,'l');
model=changeRxnBounds(model,'EX_mnl(e)',-0.148986,'l');
model=changeRxnBounds(model,'EX_fuc_L(e)',-0.148986,'l');
model=changeRxnBounds(model,'EX_glcn(e)',-0.148986,'l');
model=changeRxnBounds(model,'EX_rmn(e)',-0.148986,'l');
model=changeRxnBounds(model,'EX_arab_L(e)',-0.178783,'l');
model=changeRxnBounds(model,'EX_drib(e)',-0.178783,'l');
model=changeRxnBounds(model,'EX_rib_D(e)',-0.178783,'l');
model=changeRxnBounds(model,'EX_xyl_D(e)',-0.178783,'l');
model=changeRxnBounds(model,'EX_oxa(e)',-0.446957,'l');
model=changeRxnBounds(model,'EX_lcts(e)',-0.074493,'l');
model=changeRxnBounds(model,'EX_malt(e)',-0.074493,'l');
model=changeRxnBounds(model,'EX_sucr(e)',-0.074493,'l');
model=changeRxnBounds(model,'EX_melib(e)',-0.074493,'l');
model=changeRxnBounds(model,'EX_cellb(e)',-0.074493,'l');
model=changeRxnBounds(model,'EX_tre(e)',-0.074493,'l');
model=changeRxnBounds(model,'EX_strch1(e)',-0.257339,'l');
% fiber
model=changeRxnBounds(model,'EX_amylopect900(e)',-0.0000156731,'l');
model=changeRxnBounds(model,'EX_amylose300(e)',-0.0000470194,'l');
model=changeRxnBounds(model,'EX_arabinan101(e)',-0.0001662770,'l');
model=changeRxnBounds(model,'EX_arabinogal(e)',-0.0000219148,'l');
model=changeRxnBounds(model,'EX_arabinoxyl(e)',-0.0003066486,'l');
model=changeRxnBounds(model,'EX_bglc(e)',-0.0000000705,'l');
model=changeRxnBounds(model,'EX_cellul(e)',-0.0000282117,'l');
model=changeRxnBounds(model,'EX_dextran40(e)',-0.0001763229,'l');
model=changeRxnBounds(model,'EX_galmannan(e)',-0.0000141058,'l');
model=changeRxnBounds(model,'EX_glcmannan(e)',-0.0000328807,'l');
model=changeRxnBounds(model,'EX_homogal(e)',-0.0001282349,'l');
model=changeRxnBounds(model,'EX_inulin(e)',-0.0004701945,'l');
model=changeRxnBounds(model,'EX_kestopt(e)',-0.0028211667,'l');
model=changeRxnBounds(model,'EX_levan1000(e)',-0.0000141058,'l');
model=changeRxnBounds(model,'EX_lmn30(e)',-0.0004701945,'l');
model=changeRxnBounds(model,'EX_lichn(e)',-0.0000829755,'l');
model=changeRxnBounds(model,'EX_pect(e)',-0.0000333866,'l');
model=changeRxnBounds(model,'EX_pullulan1200(e)',-0.0000117549,'l');
model=changeRxnBounds(model,'EX_raffin(e)',-0.0047019445,'l');
model=changeRxnBounds(model,'EX_rhamnogalurI(e)',-0.0000144923,'l');
model=changeRxnBounds(model,'EX_rhamnogalurII(e)',-0.0002669874,'l');
model=changeRxnBounds(model,'EX_starch1200(e)',-0.0000117549,'l');
model=changeRxnBounds(model,'EX_xylan(e)',-0.0000320587,'l');
model=changeRxnBounds(model,'EX_xyluglc(e)',-0.0000131462,'l');
% fat
model=changeRxnBounds(model,{'EX_arachd(e)'},-0.003328,'l');
model=changeRxnBounds(model,{'EX_chsterol(e)'},-0.004958,'l');
model=changeRxnBounds(model,{'EX_glyc(e)'},-1.799655,'l');
model=changeRxnBounds(model,{'EX_hdca(e)'},-0.396371,'l');
model=changeRxnBounds(model,{'EX_hdcea(e)'},-0.036517,'l');
model=changeRxnBounds(model,{'EX_lnlc(e)'},-0.359109,'l');
model=changeRxnBounds(model,{'EX_lnlnca(e)'},-0.017565,'l');
model=changeRxnBounds(model,{'EX_lnlncg(e)'},-0.017565,'l');
model=changeRxnBounds(model,{'EX_ocdca(e)'},-0.169283,'l');
model=changeRxnBounds(model,{'EX_ocdcea(e)'},-0.681445,'l');
model=changeRxnBounds(model,{'EX_octa(e)'},-0.012943,'l');
model=changeRxnBounds(model,{'EX_ttdca(e)'},-0.068676,'l');
% protein
model=changeRxnBounds(model,{'EX_ala_L(e)';'EX_ser_L(e)';'EX_cys_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_arg_L(e)';'EX_ile_L(e)';'EX_leu_L(e)';'EX_lys_L(e)';'EX_his_L(e)'},-0.15,'l');
model=changeRxnBounds(model,{'EX_asn_L(e)';'EX_asp_L(e)';'EX_thr_L(e)'},-0.225,'l');
model=changeRxnBounds(model,{'EX_glu_L(e)';'EX_met_L(e)';'EX_gln_L(e)';'EX_pro_L(e)';'EX_val_L(e)'},-0.18,'l');
model=changeRxnBounds(model,{'EX_phe_L(e)';'EX_tyr_L(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_gly(e)'},-0.45,'l');
model=changeRxnBounds(model,{'EX_trp_L(e)'},-0.08182,'l');

% minerals and vitamins
% bacterial requirements pooled
model=changeRxnBounds(model,{'EX_thm(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ribflv(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_nac(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_btn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pnto_R(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pydxn(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pydx(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pydam(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fol(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_adpcbl(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pheme(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cbl1(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_chol(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_h2s(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_so4(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_spmd(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_na1(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_ca2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cobalt2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cl(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_k(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fe2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fe3(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_fe3dcit(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mg2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mn2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mobd(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_cu2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_zn2(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_sel(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_pime(e)'},-1,'l');
model=changeRxnBounds(model,{'EX_mqn8(e)'},-1,'l');
model=changeRxnBounds(model,'EX_12dgr180(e)',-1,'l');
model=changeRxnBounds(model,'EX_26dap_M(e)',-1,'l');
model=changeRxnBounds(model,'EX_2obut(e)',-1,'l');
model=changeRxnBounds(model,'EX_3mop(e)',-1,'l');
model=changeRxnBounds(model,'EX_4hbz(e)',-1,'l');
model=changeRxnBounds(model,'EX_ac(e)',-1,'l');
model=changeRxnBounds(model,'EX_acgam(e)',-1,'l');
model=changeRxnBounds(model,'EX_acnam(e)',-1,'l');
model=changeRxnBounds(model,'EX_acmana(e)',-1,'l');
model=changeRxnBounds(model,'EX_ade(e)',-1,'l');
model=changeRxnBounds(model,'EX_adn(e)',-1,'l');
model=changeRxnBounds(model,'EX_ala_D(e)',-1,'l');
model=changeRxnBounds(model,'EX_arab_D(e)',-1,'l');
model=changeRxnBounds(model,'EX_cgly(e)',-1,'l');
model=changeRxnBounds(model,'EX_chor(e)',-1,'l');
model=changeRxnBounds(model,'EX_cit(e)',-1,'l');
model=changeRxnBounds(model,'EX_amet(e)',-1,'l');
model=changeRxnBounds(model,'EX_csn(e)',-1,'l');
model=changeRxnBounds(model,'EX_cytd(e)',-1,'l');
model=changeRxnBounds(model,'EX_dad_2(e)',-1,'l');
model=changeRxnBounds(model,'EX_dcyt(e)',-1,'l');
model=changeRxnBounds(model,'EX_ddca(e)',-1,'l');
model=changeRxnBounds(model,'EX_dgsn(e)',-1,'l');
model=changeRxnBounds(model,'EX_for(e)',-1,'l');
model=changeRxnBounds(model,'EX_fum(e)',-1,'l');
model=changeRxnBounds(model,'EX_glyc3p(e)',-1,'l');
model=changeRxnBounds(model,'EX_gsn(e)',-1,'l');
model=changeRxnBounds(model,'EX_gthox(e)',-1,'l');
model=changeRxnBounds(model,'EX_gua(e)',-1,'l');
model=changeRxnBounds(model,'EX_h(e)',-1,'l');
model=changeRxnBounds(model,'EX_hxan(e)',-1,'l');
model=changeRxnBounds(model,'EX_indole(e)',-1,'l');
model=changeRxnBounds(model,'EX_lanost(e)',-1,'l');
model=changeRxnBounds(model,'EX_metsox_S_L(e)',-1,'l');
model=changeRxnBounds(model,'EX_mqn7(e)',-1,'l');
model=changeRxnBounds(model,'EX_ncam(e)',-1,'l');
model=changeRxnBounds(model,'EX_nmn(e)',-1,'l');
model=changeRxnBounds(model,'EX_orn(e)',-1,'l');
model=changeRxnBounds(model,'EX_pheme(e)',-1,'l');
model=changeRxnBounds(model,'EX_pi(e)',-1,'l');
model=changeRxnBounds(model,'EX_ptrc(e)',-1,'l');
model=changeRxnBounds(model,'EX_pydx5p(e)',-1,'l');
model=changeRxnBounds(model,'EX_q8(e)',-1,'l');
model=changeRxnBounds(model,'EX_sheme(e)',-1,'l');
model=changeRxnBounds(model,'EX_spmd(e)',-1,'l');
model=changeRxnBounds(model,'EX_thm(e)',-1,'l');
model=changeRxnBounds(model,'EX_thymd(e)',-1,'l');
model=changeRxnBounds(model,'EX_ura(e)',-1,'l');
model=changeRxnBounds(model,'EX_uri(e)',-1,'l');
model=changeRxnBounds(model,'EX_xan(e)',-1,'l');
model=changeRxnBounds(model,'EX_2dmmq8(e)',-1,'l');
model=changeRxnBounds(model,'EX_4abz(e)',-1,'l');
model=changeRxnBounds(model,'EX_amp(e)',-1,'l');
model=changeRxnBounds(model,'EX_fald(e)',-1,'l');
model=changeRxnBounds(model,'EX_gam(e)',-1,'l');
model=changeRxnBounds(model,'EX_glu_D(e)',-1,'l');
model=changeRxnBounds(model,'EX_gthrd(e)',-1,'l');
model=changeRxnBounds(model,'EX_h2(e)',-1,'l');
model=changeRxnBounds(model,'EX_no2(e)',-1,'l');
% only for methanogens
model=changeRxnBounds(model,'EX_meoh(e)',-10,'l');
% other compounds
model=changeRxnBounds(model,'EX_h2o(e)',-10,'l');

modelOut=model;

end