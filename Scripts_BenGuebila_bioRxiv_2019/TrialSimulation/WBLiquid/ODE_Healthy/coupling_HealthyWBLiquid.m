function coupling_HealthyWBLiquid(simType,offTarget)
    %cd(pathOs('p\Projects\Whole_body_model\GIM\TrialSimulation'));
    load simResultsTrials.mat;
    yout = simResultsTrials{13};%umol 
    diseaseState = 'Healthy';
    trialCondition = 'WBLiquid';
    simLauncher;
end