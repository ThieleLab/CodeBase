%
%script to generate time points of ODE simulations
% and batch run coupling algorithm
% Marouen BEN GUEBILA 10/2015
%%
trialConditions = {'IVITT' 'IVITT' 'IVGTT' 'IVGTT' 'Noinf' 'Noinf' 'SCIB' 'SCIB' 'SCII' 'SCII' 'WBSolid' 'WBSolid' 'WBLiquid' 'WBLiquid'};
simResultsTrials = cell(length(trialConditions)+1,1);
simResultsTrials{length(trialConditions)+1} = [0:5:600]; 
parpool(14);
%%
parfor i=1:length(trialConditions)
    if mod(i,2)
        %IVITT
        cd(pathOs(['..' filesep 'TrialSimulation' filesep trialConditions{i} filesep 'ODE_Healthy']));
        diseaseState = 'Healthy';
        trialCondition = trialConditions{i};
        simName = ['ODEMain_' diseaseState trialCondition]
        fh = str2func(simName);
        [tout, yout] = fh(); 
        yout = yout(11:length(tout),:);%umol
        tout = tout(11:length(tout))-1000;
        simResultsTrials{i} = yout;
    else
        cd(pathOs(['..' filesep 'TrialSimulation' filesep trialConditions{i} filesep 'ODE_T1D']));
        diseaseState = 'T1D';
        trialCondition = trialConditions{i};
        simName = ['ODEMain_' diseaseState trialCondition]
        fh = str2func(simName);
        [tout, yout] = fh(); 
        yout = yout(11:length(tout),:);%umol
        tout = tout(11:length(tout))-1000;
        simResultsTrials{i} = yout;
    end
end
cd(['..' filesep '..' filesep '..' filesep 'data'])
save('simResultsTrials.mat','simResultsTrials');