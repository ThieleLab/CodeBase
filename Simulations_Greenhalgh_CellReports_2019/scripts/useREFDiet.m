function model=useREFDiet(model)
% Implements the REF diet (fiber-free reference diet) as constraints in the
% in silico models.

model=changeRxnBounds(model,model.rxns(strmatch('EX_',model.rxns)),0,'l');
%% amino acids
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
model=changeRxnBounds(model,'EX_ala_L[u]',-1,'l');
model=changeRxnBounds(model,'EX_glc_D[u]',-10,'l');
%% other
model=changeRxnBounds(model,'EX_pyr[u]',-1,'l');
model=changeRxnBounds(model,'EX_lnlc[u]',-1,'l');
model=changeRxnBounds(model,'EX_lipoate[u]',-1,'l');

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
