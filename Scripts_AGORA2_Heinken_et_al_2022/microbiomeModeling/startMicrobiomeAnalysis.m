% start microbiome analysis

objectiveList={
    'EX_sn38[fe]'   'Diet_EX_sn38g[d]'
    'EX_r406[fe]'   'Diet_EX_r788[d]'
    'EX_5fura[fe]'   'Diet_EX_fcsn[d]'
    'EX_dh5fura[fe]'   'Diet_EX_fcsn[d]'
    'EX_dh5fura[fe]'   'Diet_EX_5fura[d]'
    'EX_dfduri[fe]' 'Diet_EX_dfdcytd[d]'
    'EX_dihydro_digoxin[fe]'    'Diet_EX_digoxin[d]'
    'EX_ac5asa[fe]' 'Diet_EX_5asa[d]'
    'EX_5asa[fe]'   'Diet_EX_bzd[d]'
    'EX_ac5asa[fe]'   'Diet_EX_bzd[d]'
    'EX_nchlphncl[fe]'  'Diet_EX_chlphncl[d]'
    'EX_bvu[fe]'    'Diet_EX_srv[d]'
    'EX_dopa[fe]'   'Diet_EX_34dhphe[d]'
    'EX_mtym[fe]'   'Diet_EX_34dhphe[d]'
    'EX_pcresol[fe]'    'Diet_EX_4hphac[d]'
    'EX_cholate[fe]'    'Diet_EX_tchola[d]'
    };

%% Japanese diet
modelFolder=[rootDir filesep 'Modeling_CRC' filesep 'JapaneseDiet'];
solutionFolder=[rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices_JD' filesep];

%% Run the computation
[objectives,shadowPrices]=analyseObjectiveShadowPrices(modelFolder,objectiveList,'SPDef','Nonzero','numWorkers',numWorkers,'resultsFolder',solutionFolder);
objectives=cell2table(objectives');
writetable(objectives,[solutionFolder 'AGORA2_CRC_Objectives_JD'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

shadowPrices=cell2table(shadowPrices);
writetable(shadowPrices,[solutionFolder 'AGORA2_CRC_ShadowPrices_JD'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

%% Average European diet
modelFolder=[rootDir filesep 'Modeling_CRC' filesep 'AverageEuropeanDiet'];
solutionFolder=[rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices_AED' filesep];

%% Run the computation
[objectives,shadowPrices]=analyseObjectiveShadowPrices(modelFolder,objectiveList,'SPDef','Nonzero','numWorkers',numWorkers,'resultsFolder',solutionFolder);
objectives=cell2table(objectives');
writetable(objectives,[solutionFolder 'AGORA2_CRC_Objectives_AED'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

shadowPrices=cell2table(shadowPrices);
writetable(shadowPrices,[solutionFolder 'AGORA2_CRC_ShadowPrices_AED'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
