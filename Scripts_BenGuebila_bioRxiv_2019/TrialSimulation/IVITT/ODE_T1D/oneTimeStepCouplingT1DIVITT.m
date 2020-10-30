function [FBA,harvey] = oneTimeStepCouplingT1DIVITT(i,harvey,yout,tout,diseaseState,trialCondition,simType,offTarget,prevTSFD)
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
    end
    %%
    %Change objective function
    harvey = changeObjective(harvey,'Whole_body_objective_rxn');
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