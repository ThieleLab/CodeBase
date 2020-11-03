%%
%fastFVA
parpool(32);
changeCobraSolver('tomlab_cplex');
cpxControl = CPLEXParamSet_Harvey;
load augHarvey;
harvey = modelOrganAllCoupled;
%access save folder
cd(['..' filesep '..' filesep '..' filesep 'data' filesep 'FVA'])
%%
%subject constraints
stdCons;
%%
diseaseState = 'Unc';
[minFluxUnc,maxFluxUnc,optsolUnc,retUnc,fbasolUnc,fvaminUnc,fvamaxUnc,statussolmin,statussolmax] =...
    fastFVA(harvey,90,'max',[],[],'A',cpxControl);
save(['minFlux' diseaseState '.mat'],'minFluxUnc','-v7.3');
save(['maxFlux' diseaseState '.mat'],'maxFluxUnc','-v7.3');
save(['fvamin' diseaseState '.mat'],'fvaminUnc','-v7.3');
save(['fvamax' diseaseState '.mat'],'fvamaxUnc','-v7.3');