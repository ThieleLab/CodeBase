function coupling_HealthyWBSolid(simType,offTarget)
    %cd(pathOs('p\Projects\Whole_body_model\GIM\TrialSimulation'));
    load simResultsTrials.mat;
    yout = simResultsTrials{11};%umol 
    diseaseState = 'Healthy';
    trialCondition = 'WBSolid';
    simLauncher;
end