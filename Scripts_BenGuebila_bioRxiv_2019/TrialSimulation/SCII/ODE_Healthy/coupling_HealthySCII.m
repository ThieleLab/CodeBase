function coupling_HealthySCII(simType,offTarget)
    %cd(pathOs('p\Projects\Whole_body_model\GIM\TrialSimulation'));
    load simResultsTrials.mat;
    yout = simResultsTrials{9};%umol 
    diseaseState = 'Healthy';
    trialCondition = 'SCII';
    simLauncher;
end