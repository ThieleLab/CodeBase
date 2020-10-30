function coupling_T1DSCII(simType,offTarget)
    load simResultsTrials.mat;
    yout = simResultsTrials{10};%umol 
    diseaseState = 'T1D';
    trialCondition = 'SCII';
    simLauncher;
end