function coupling_T1DWBSolid(simType,offTarget)
    load simResultsTrials.mat;
    yout = simResultsTrials{12};%umol 
    diseaseState = 'T1D';
    trialCondition = 'WBSolid';
    simLauncher;
end