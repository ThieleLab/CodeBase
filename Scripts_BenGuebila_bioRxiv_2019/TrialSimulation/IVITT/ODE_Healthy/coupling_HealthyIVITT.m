function coupling_HealthyIVITT(simType,offTarget)
    %cd(pathOs('p\Projects\Whole_body_model\GIM\TrialSimulation'));
    load simResultsTrials.mat;
    yout = simResultsTrials{1};%umol 
    diseaseState = 'Healthy';
    trialCondition = 'IVITT';
    simLauncher;
end