%%
% Coupling Harvey and GIM
% Glucose insulin Glucagon regulatory and metabolic system
% BEN GUEBILA Marouen 05/08/2015
%%
clear harvey modelOrganAllCoupled
load('augHarvey')
harvey_model = modelOrganAllCoupled;
sprintf('///////////////////////New simulation////////////////////////')
% GIM
%Noinf
tout = simResultsTrials{15};
%%
tic
% Coupling loop;
harvey_model = changeObjective(harvey_model,'Whole_body_objective_rxn');
global restart;
restart
if ~isempty(restart) && isequal(restart,1)
    load(['fluxMean' diseaseState trialCondition simType offTarget 'P' num2str(patient) '.mat']);
    load(['errVec' diseaseState trialCondition simType offTarget 'P' num2str(patient) '.mat']);
    load(['objVec' diseaseState trialCondition simType offTarget 'P' num2str(patient) '.mat']);
    load(['obj' diseaseState trialCondition simType  offTarget 'P' num2str(patient) '.mat']);
    indNnz = find(cellfun(@isempty,objVec));
    start=indNnz(1);
    FBA.full=objVec{start-1};
else
    start=1;
    FBA.full=[];
    errVec     = zeros(length(tout)-1,1);
    fluxMean   = zeros(length(tout)-1,1);
    objVec     = cell(length(tout)-1,1);
    obj = objVec;
end
if or(isequal(offTarget,'GenExIns1h'),isequal(offTarget,'GenExIns1hDMLNAA'))
    load timeInsulin5var;
    tout = timeInsulin5var;
end
for i = 1:length(tout)-1
    OneStepFunString = ['oneTimeStepCoupling' diseaseState trialCondition];
    OneStepFunName = str2func(OneStepFunString)
    if isequal(simType,'CRONICS')
        [FBA,harvey] = OneStepFunName(i,harvey_model,yout,tout,diseaseState,trialCondition,simType,offTarget,FBA.full);
    else
        [FBA,harvey] = OneStepFunName(i,harvey_model,yout,tout,diseaseState,trialCondition,simType,offTarget);
    end
    errVec(i) = FBA.stat;
    objVec{i} = FBA.full;%changed 
    fluxMean(i) = mean(abs(FBA.full));
    obj{i} = FBA.obj;
    if mod(i,10)==0
        if exist('patient')
            cd(['patients' filesep 'results'])
            save(['fluxMean' diseaseState trialCondition simType offTarget 'P' num2str(patient) '.mat'],'fluxMean')
            save(['errVec' diseaseState trialCondition simType offTarget 'P' num2str(patient) '.mat'],'errVec')
            save(['objVec' diseaseState trialCondition simType offTarget 'P' num2str(patient) '.mat'],'objVec')
            save(['obj' diseaseState trialCondition simType  offTarget 'P' num2str(patient) '.mat'],'obj')
            cd ../..
        else
            cd(['results_' diseaseState '_' trialCondition])
            save(['fluxMean' diseaseState trialCondition simType offTarget '.mat'],'fluxMean')
            save(['errVec' diseaseState trialCondition simType offTarget '.mat'],'errVec')
            save(['objVec' diseaseState trialCondition simType offTarget '.mat'],'objVec')
            save(['obj' diseaseState trialCondition simType  offTarget '.mat'],'obj')
            cd ..
        end
    end
end

if exist('patient')
    cd(['patients' filesep 'results'])
    save(['fluxMean' diseaseState trialCondition simType offTarget 'P' num2str(patient) '.mat'],'fluxMean')
    save(['errVec' diseaseState trialCondition simType offTarget 'P' num2str(patient) '.mat'],'errVec')
    save(['objVec' diseaseState trialCondition simType offTarget 'P' num2str(patient) '.mat'],'objVec')
    save(['obj' diseaseState trialCondition simType  offTarget 'P' num2str(patient) '.mat'],'obj')
else
    cd(['results_' diseaseState '_' trialCondition])
    save(['fluxMean' diseaseState trialCondition simType offTarget '.mat'],'fluxMean')
    save(['errVec' diseaseState trialCondition simType offTarget '.mat'],'errVec')
    save(['objVec' diseaseState trialCondition simType offTarget '.mat'],'objVec')
    save(['obj' diseaseState trialCondition simType  offTarget '.mat'],'obj')
end
toc