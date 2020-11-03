function coupling_T1DNoinf(simType,offTarget)
    load simResultsTrials.mat;
    yout = simResultsTrials{6};%umol 
    diseaseState = 'T1D';
    trialCondition = 'Noinf';
    simLauncher;
end