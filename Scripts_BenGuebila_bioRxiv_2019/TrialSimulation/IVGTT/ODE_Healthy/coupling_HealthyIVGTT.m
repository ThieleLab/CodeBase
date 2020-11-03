function coupling_HealthyIVGTT(simType,offTarget)
    load simResultsTrials.mat;
    yout = simResultsTrials{3};%umol 
    diseaseState = 'Healthy';
    trialCondition = 'IVGTT';
    simLauncher;
end
