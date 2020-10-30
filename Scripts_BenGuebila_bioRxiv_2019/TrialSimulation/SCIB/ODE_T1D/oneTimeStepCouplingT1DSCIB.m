function [FBA,harvey] = oneTimeStepCouplingT1DSCIB(i,harvey,yout,tout,diseaseState,trialCondition,simType,offTarget,prevTSFD)
    i
    %%
    %convert constraints to  /5mn
    indu = find(harvey.ub ~= 1000000 & harvey.ub~= 0 & harvey.ub~= -1000000);
    indl = find(harvey.lb ~= 1000000 & harvey.lb~= 0 & harvey.lb~= -1000000);
    harvey.lb(indl) = harvey.lb(indl)/24/60*5;
    harvey.ub(indu) = harvey.ub(indu)/24/60*5;
    %%
    if isequal(offTarget,'GenEx')
        updateConstraintsEx;
    elseif isequal(offTarget,'GenExIns1h')
        updateConstraintsEx;
        updateConstraintsIns1h;
    elseif or(isequal(offTarget,'GenExInsTot'),isequal(offTarget,'GenExInsTotDMLNAA'))
        updateConstraintsEx;
        updateConstraintsInsTot;
    end
    %%
    %Change objective function
    if or(isequal(offTarget,'GenExIns1hDMLNAA'),isequal(offTarget,'GenExInsTotDMLNAA'))
        [harvey,rxnNames] = addDemandReaction(harvey,{'phe_L[bc]','val_L[bc]','met_L[bc]','ile_L[bc]','leu_L[bc]','tyr_L[bc]'...
        'his_L[bc]','trp_L[bc]','lys_L[bc]'});
        harvey.c   = zeros(length(harvey.c),1);
        harvey.c(end-8:end) = ones(9,1);%LNAA demand reaction
    else
        harvey = changeObjective(harvey,'Whole_body_objective_rxn');
    end
    %%
    %Parameters
    Time = 1000+tout(i);
    y = yout(i,:);
    ODERHSFunction(1000+tout(i), y);
    [y switchUpdate] = PerformSwitches(Time, y);
    paramString = ['ParametersGIM_' trialCondition '_' diseaseState];
    paramStringFunName = str2func(paramString);
    paramStringFunName();
    %%
    oneTimeStepConstraints;
    %%
    %Optimize model
    optimModelT1D;
end