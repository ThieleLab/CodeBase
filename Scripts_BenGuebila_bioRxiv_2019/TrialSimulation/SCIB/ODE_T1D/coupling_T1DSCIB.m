function coupling_T1DSCIB(simType,offTarget,patient)
    if nargin>2
        load simResultsTrials.mat;%for time vector
        load simRecord31.mat;
        yout = simRecord{patient};%umol 
    else
        %cd(pathOs('p\Projects\Whole_body_model\GIM\TrialSimulation'));
        load simResultsTrials.mat;
        yout = simResultsTrials{8};%umol 
    end
    diseaseState = 'T1D';
    trialCondition = 'SCIB';
    simLauncher;
end