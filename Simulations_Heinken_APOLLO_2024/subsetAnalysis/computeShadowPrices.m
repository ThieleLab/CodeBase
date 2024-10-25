
% compute production of key metabolites
% modPath=[rootDir filesep 'MicrobiomeModeling' filesep 'Microbiome_models_AED'];
modPath=[pwd filesep 'Microbiome_models_AED'];

numWorkers = 12;

% some potentially interesting metabolites
metList={'EX_ac[fe]'
'EX_ppa[fe]'
'EX_but[fe]'
'EX_isobut[fe]'
'EX_isoval[fe]'
'EX_lac_D[fe]'
'EX_lac_L[fe]'
'EX_for[fe]'
'EX_etoh[fe]'
'EX_h2s[fe]'
'EX_tma[fe]'
'EX_phe_L[fe]'
'EX_tyr_L[fe]'
'EX_trp_L[fe]'
'EX_dopa[fe]'
'EX_taur[fe]'
'EX_pcresol[fe]'
'EX_indole[fe]'
'EX_4abut[fe]'
'EX_leu_L[fe]'
'EX_ile_L[fe]'
'EX_val_L[fe]'
'EX_cholate[fe]'
'EX_dchac[fe]'
'EX_HC02191[fe]'
'EX_12dhchol[fe]'
'EX_7ocholate[fe]'};

% spPath = [rootDir filesep 'MicrobiomeModeling' filesep 'ShadowPrices'];
spPath = [pwd filesep 'ShadowPrices'];

[objectives,shadowPrices]=analyseObjectiveShadowPrices(modPath, metList, 'resultsFolder', spPath, 'numWorkers', numWorkers);

% writetable(cell2table(objectives),[rootDir filesep 'MicrobiomeModeling' filesep 'Objectives_AED'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
writetable(cell2table(objectives),[pwd filesep 'Objectives_AED'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

shadowPrices(1,:)=strrep(shadowPrices(1,:),'.mat','');
% writetable(cell2table(shadowPrices),[rootDir filesep 'MicrobiomeModeling' filesep 'ShadowPrices_AED'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
writetable(cell2table(shadowPrices),[pwd filesep 'ShadowPrices_AED'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
