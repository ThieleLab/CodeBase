function coupling_T1DIVITT(simType,offTarget)
    load simResultsTrials.mat;
    yout = simResultsTrials{2};%umol 
    diseaseState = 'T1D';
    trialCondition = 'IVITT';
    simLauncher;
end
