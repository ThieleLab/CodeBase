function model=useHFDiet(model)
% Implements the HF diet (high-fiber diet) as constraints in the
% in silico models.

% casein & bactopeptone
model=changeRxnBounds(model,model.rxns(strmatch('EX_',model.rxns)),0,'l');
model=changeRxnBounds(model,'EX_gly[u]',-1,'l');
model=changeRxnBounds(model,'EX_arg_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_cys_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_gln_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_his_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_ile_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_leu_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_lys_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_met_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_phe_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_ser_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_thr_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_trp_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_tyr_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_val_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_glu_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_pro_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_ala_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_asn_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_asp_L[u]',-1,'l');

% ox bile=bile acids
model=changeRxnBounds(model,'EX_cholate[u]',-1,'l');
model=changeRxnBounds(model,'EX_C02528[u]',-1,'l');
model=changeRxnBounds(model,'EX_tchola[u]',-1,'l');
model=changeRxnBounds(model,'EX_gchola[u]',-1,'l');
model=changeRxnBounds(model,'EX_dgchol[u]',-1,'l');
model=changeRxnBounds(model,'EX_tdchola[u]',-1,'l');

% fibers
model=changeRxnBounds(model,'EX_pect[u]',-0.000094,'l');
model=changeRxnBounds(model,'EX_xylan[u]',-0.00094,'l');
model=changeRxnBounds(model,'EX_arabinogal[u]',-0.00094,'l');
model=changeRxnBounds(model,'EX_amylopect900[u]',-0.00094,'l');
model=changeRxnBounds(model,'EX_strch1[u]',-0.094,'l');

%% ions and vitamins: no exact composition given in the experimental medium
% I assume the ones required by the reconstructions
% ions
model=changeRxnBounds(model,'EX_ca2[u]',-1,'l');
model=changeRxnBounds(model,'EX_cl[u]',-1,'l');
model=changeRxnBounds(model,'EX_so4[u]',-1,'l');
model=changeRxnBounds(model,'EX_h2s[u]',-1,'l');
model=changeRxnBounds(model,'EX_cobalt2[u]',-1,'l');
model=changeRxnBounds(model,'EX_cu2[u]',-1,'l');
model=changeRxnBounds(model,'EX_fe2[u]',-1,'l');
model=changeRxnBounds(model,'EX_fe3[u]',-1,'l');
model=changeRxnBounds(model,'EX_k[u]',-1,'l');
model=changeRxnBounds(model,'EX_mg2[u]',-1,'l');
model=changeRxnBounds(model,'EX_mn2[u]',-1,'l');
model=changeRxnBounds(model,'EX_zn2[u]',-1,'l');
model=changeRxnBounds(model,'EX_pi[u]',-10,'l');
model=changeRxnBounds(model,'EX_h2o[u]',-10,'l');

%% vitamins
model=changeRxnBounds(model,'EX_fol[u]',-1,'l');
model=changeRxnBounds(model,'EX_inost[u]',-1,'l');
model=changeRxnBounds(model,'EX_nac[u]',-1,'l');
model=changeRxnBounds(model,'EX_ncam[u]',-1,'l');
model=changeRxnBounds(model,'EX_pnto_R[u]',-1,'l');
model=changeRxnBounds(model,'EX_pydx[u]',-1,'l');
model=changeRxnBounds(model,'EX_pydxn[u]',-1,'l');
model=changeRxnBounds(model,'EX_ribflv[u]',-1,'l');
model=changeRxnBounds(model,'EX_sheme[u]',-1,'l');
model=changeRxnBounds(model,'EX_thm[u]',-1,'l');

% supplemented with hemin and vitamin K
model=changeRxnBounds(model,'EX_pheme[u]',-1,'l');
model=changeRxnBounds(model,'EX_mqn7[u]',-1,'l');
model=changeRxnBounds(model,'EX_mqn8[u]',-1,'l');
model=changeRxnBounds(model,'EX_q8[u]',-1,'l');

end
