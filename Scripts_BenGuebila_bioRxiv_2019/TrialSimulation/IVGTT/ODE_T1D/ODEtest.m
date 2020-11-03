cd(pathOs('P:\Projects\Whole_body_model\GIM\TrialSimulation'));    
load simResultsTrials.mat;
yout = simResultsTrials{4};%umol 
tout = simResultsTrials{15};
diseaseState = 'T1D';
trialCondition = 'IVGTT';
cd(pathOs(['P:\Projects\Whole_body_model\GIM\TrialSimulation\' trialCondition '\ODE_' diseaseState]))
%Parameters
Yout=[];
for i=1:length(tout)
    i
    Time = 1000+tout(i);
    y = yout(i,:);
    ODERHSFunction(1000+tout(i), y);
    [y switchUpdate] = PerformSwitches(Time, y);
    paramString = ['ParametersGIM_' trialCondition '_' diseaseState];
    paramStringFunName = str2func(paramString);
    paramStringFunName();
    Yout = [Yout y(1)];
end
plot(tout,Yout,tout,yout(:,1));
legend('ODEtest','OriginalSim')