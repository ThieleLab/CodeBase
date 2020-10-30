function coupling_T1DWBLiquid(simType,offTarget)
    load simResultsTrials.mat;
    yout = simResultsTrials{14};%umol 
    diseaseState = 'T1D';
    trialCondition = 'WBLiquid';
    simLauncher;
end