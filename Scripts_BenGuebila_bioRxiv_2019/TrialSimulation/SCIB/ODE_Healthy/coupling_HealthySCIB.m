function coupling_HealthySCIB(simType,offTarget)
    %cd(pathOs('p\Projects\Whole_body_model\GIM\TrialSimulation'));
    load simResultsTrials.mat;
    yout = simResultsTrials{7};%umol 
    diseaseState = 'Healthy';
    trialCondition = 'SCIB';
    simLauncher;
end