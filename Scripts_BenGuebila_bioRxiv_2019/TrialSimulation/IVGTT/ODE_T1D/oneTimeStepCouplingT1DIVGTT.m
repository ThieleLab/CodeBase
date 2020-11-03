function [FBA,harvey] = oneTimeStepCouplingT1DIVGTT(i,harvey,yout,tout,diseaseState,trialCondition,simType,offTarget,prevTSFD)
    i
    %%
    %convert constraints to  /5mn
    indu = find(harvey.ub ~= 1000000 & harvey.ub~= 0 & harvey.ub~= -1000000);
    indl = find(harvey.lb ~= 1000000 & harvey.lb~= 0 & harvey.lb~= -1000000);
    harvey.lb(indl) = harvey.lb(indl)/24/60*5;
    harvey.ub(indu) = harvey.ub(indu)/24/60*5;
    %%
    if isequal(offTarget,'GenEx') || isequal(offTarget,'GenExDMATP')
        updateConstraintsEx;
    end
    %%
    %Change objective function
    if isequal(offTarget,'GenExDMATP')
        harvey.c = zeros(length(harvey.c),1);
        indAtpDM = find(cellfun(@(x) ~isempty(strfind(x,'DM_atp_c_')),harvey.rxns));
        harvey.c(indAtpDM) = 1;
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