function coupling_T1DIVGTT(simType,offTarget)
    load simResultsTrials.mat;
    yout = simResultsTrials{4};%umol 
    diseaseState = 'T1D';
    trialCondition = 'IVGTT';
    simLauncher;
end