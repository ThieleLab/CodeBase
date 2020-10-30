%%
%fastFVA
parpool(16);
changeCobraSolver('tomlab_cplex');
cpxControl = CPLEXParamSet_Harvey;
load augHarvey;
harvey = modelOrganAllCoupled;
%%
%Reversible model
modelType      = 'rev';
trialCondition = 'Noinf';
diseaseState   = 'T1D';
%access save folder
cd(['..' filesep '..' filesep '..' filesep 'data' filesep 'FVA'])
load simResultsTrials.mat;
%subjects constraint
[harvey,irrevMatch] = convertHarveyIrrev(harvey,modelType,trialCondition,diseaseState,simResultsTrials,cpxControl);
%%
[minFluxH,maxFluxH,optsolUnc,retUnc,fbasolUnc,fvaminH,fvamaxH,statussolmin,statussolmax] =...
    fastFVA(harvey,90,'max',[],[],'A',cpxControl);
%access save folder
cd(['..' filesep '..' filesep '..' filesep 'data' filesep 'FVA'])
save(['minFlux' diseaseState modelType '.mat'],'minFluxH','-v7.3');
save(['maxFlux' diseaseState modelType '.mat'],'maxFluxH','-v7.3');
save(['fvamin' diseaseState modelType '.mat'],'fvaminH','-v7.3');
save(['fvamax' diseaseState modelType '.mat'],'fvamaxH','-v7.3');
%%
%seperatly save reactions needed for insulin constraints
indRxn=findRxnIDs(harvey,{'Liver_PGI','Liver_TALA','Liver_FBP','Liver_PFK',...
    'Liver_FBA','Liver_TPI','Liver_PGK','Liver_PGM','Liver_ENO'});
fvamaxTI = fvamaxH(indRxn,:);
fvaminTI = fvaminH(indRxn,:);
save('fvamaxTI','fvamaxTI');
save('fvaminTI','fvaminTI');
%%
%Irreversible Model
%subject constraints and sets to irreversible
modelType      = 'Irrev';
harvey = modelOrganAllCoupled;
cd(['..' filesep '..' filesep 'common' filesep 'FVA' filesep 't1d'])
[harvey,irrevMatch] = convertHarveyIrrev(harvey,modelType,trialCondition,diseaseState,simResultsTrials,cpxControl);
%%
[minFluxH,maxFluxH,optsolUnc,retUnc,fbasolUnc,fvaminH,fvamaxH,statussolmin,statussolmax] =...
    fastFVA(harvey,90,'max',[],[],'A',cpxControl);
%access save folder
cd(['..' filesep '..' filesep '..' filesep 'data' filesep 'FVA'])
save(['minFlux' diseaseState modelType '.mat'],'minFluxH','-v7.3');
save(['maxFlux' diseaseState modelType '.mat'],'maxFluxH','-v7.3');
save(['fvamin' diseaseState modelType '.mat'],'fvaminH','-v7.3');
save(['fvamax' diseaseState modelType '.mat'],'fvamaxH','-v7.3');
if isequal(modelType,'Irrev')
    save(['irrevMatch' diseaseState '.mat'],'irrevMatch','-v7.3');
end