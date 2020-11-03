function coupling_HealthyNoinf(simType,offTarget)
    load simResultsTrials.mat;
    yout = simResultsTrials{5};%umol 
    diseaseState = 'Healthy';
    trialCondition = 'Noinf';
    simLauncher;
end