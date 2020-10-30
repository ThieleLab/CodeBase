function [FBA,harvey] = oneTimeStepCouplingHealthyNoinf(i,harvey,yout,tout,diseaseState,trialCondition,simType,offTarget,prevTSFD)
    i
    %%
    %allow some organs ot secrete glucose(diffusion back)
%     oldNumVec=find(cellfun(@(x) ~isempty(strmatch(x,{...
%     'Muscle_EX_glc(e)_[bc]'...
%     'Liver_EX_glc(e)_[bc]'...
%     'Pancreas_EX_glc(e)_[bc]'...
%     'Kidney_EX_glc(e)_[bc]'...
%     'Liver_sink_glygn2(c)'})),harvey.rxns));%
%     harvey.lb(oldNumVec)=-1e+6;
%     harvey.ub(oldNumVec)=+1e+6;
%     harvey = changeRxnBounds(harvey,'EX_h2o[sw]',0,'l');
%     harvey = changeRxnBounds(harvey,'EX_h2o[sw]',+1e+6,'u');
    %%
    %convert constraints to  /5mn
    indu = find(harvey.ub ~= 1000000 & harvey.ub~= 0 & harvey.ub~= -1000000);
    indl = find(harvey.lb ~= 1000000 & harvey.lb~= 0 & harvey.lb~= -1000000);
    harvey.lb(indl) = harvey.lb(indl)/24/60*5;
    harvey.ub(indu) = harvey.ub(indu)/24/60*5;
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
    optimModel;
end