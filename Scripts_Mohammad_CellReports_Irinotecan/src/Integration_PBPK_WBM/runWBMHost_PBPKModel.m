function [stat, Flux_t, solution1,solution,t] = runWBMHost_PBPKModel(str, PBPKModelFoldername, personalizedModelFoldername)
% Simulates a patient-specific whole-body distribution model integrated with a PBPK model.
%
% Loads patient-specific models and integrates them with a physiologically based pharmacokinetic model to simulate drug distribution and metabolism.
%
% USAGE:
%
%   [stat, Flux_t, solution1, solution, t] = runWBMHost_PBPKModel(str, PBPKModelFoldername, personalizedModelFoldername)
%
% INPUT:
%    patientID:                Identifier for the patient.
%    str:                      A string used for naming conventions.
%    PBPKModelFoldername:      Folder name containing PBPK model results.
%    personalizedModelFoldername: Folder name containing personalized model results.
%    Physiological_param:      Struct containing physiological parameters for the patient.
%    Dosage_inf:               Struct containing dosage information.
%    Drug_param:               Struct containing drug parameters.
%
% OUTPUT:
%    stat:          Status of the FBA solutions.
%    Flux_t:        Time series data for fluxes.
%    solution1:     Solution of the model with zero dosage.
%    solution:      Reference solution without drug metabolites.
%    t:             Time points for the simulation.
%
% DETAILED DESCRIPTION:
% The function performs the following steps:
% 1. Loads the patient-specific whole-body distribution model.
% 2. Initializes the model and sets the objective reaction bounds.
% 3. Identifies and sets reaction IDs and bounds for specific metabolites.
% 4. Rounds model parameters to a specified precision.
% 5. Loads PBPK simulation results.
% 6. Ensures the model has no objective except the one set for optimization.
% 7. Applies PBPK constraints to the model with zero dosage and concentration.
% 8. Optimizes the model with zero dosage and concentration.
% 9. Creates a reference model without drug metabolites and optimizes it.
% 10. Converts concentration values to mmol.
% 11. Performs flux balance analysis (FBA) for each time point.
% 12. Combines FBA results.
% 13. Saves the status of the solution obtained and extracts time series data for fluxes.

% July 2024 - Faiz Khan Mohammad & Ines Thiele

% Load patient-specific Irinoitecan Whole-Body Distribution Model (WBDM)
load(strcat(['Results\',personalizedModelFoldername,'\WBM_', str, '.mat']));

% Initialize the model and set the objective reaction bounds
model.c(:) = 0;
model = changeRxnBounds(model, 'Whole_body_objective_rxn', 1, 'b');

% Find and set reaction IDs and bounds for specific metabolites
rxns_id = findRxnIDs(model, findRxnsFromMets(model, model.mets(contains(model.mets, {'cpt11[', 'sn38[', 'sn38g[', 'apc[', ...
        'npc[', 'apcald[', 'apca[','sn38cb[', 'sn38gcb[', 'cpt11cb[', 'bpca[', 'etiri[', 'tetrapeg[', 'lvlnar['}))));
model.lb(rxns_id(model.lb(rxns_id) ~= 0)) = -1000;
model.ub(rxns_id) = 1000;

% Round model parameters to a specific precision
param.rounding = 5;
model.lb = round(model.lb, param.rounding);
model.ub = round(model.ub, param.rounding);
model.S = round(model.S, param.rounding);

% Load PBPK simulation results
load(strcat(['Results\',PBPKModelFoldername,'\PBPK_simulation_', str, '.mat']));

% Ensure the model has no objective except the one set for optimization
param.rounding = 5;
param.solverName = 'ibm_cplex';
param.minNorm = 1e-5;

% Note: High performance computing might be required for optimization
% Store the original model
model1 = model;

% Initialize zero dosage and concentration
warning('off', 'all');
zero_dosage = Dosage_inf;
zero_dosage.dose_amt = 0;
zero_C = zeros(23 * Drug_param.N_metabolites, 1);

% Apply PBPK constraints to the model
model1 = PBPK_constraint_WBM(model1, Physiological_param, zero_dosage, Drug_param, zero_C, 0);
model1.osenseStr = 'min';
solution1 = optimizeWBModel(model1, param);

modelref = model;
modelref = changeRxnBounds(modelref,modelref.rxns(rxns_id),0,'b');
modelref.osenseStr = 'min';
solution = optimizeWBModel(modelref, param);

% Convert concentration values to mmol
for i = 1:size(C, 1)
    C(i, :) = C(i, :) / Drug_param.metabolite{fix((i - 1) / 23) + 1}.MW; % in mmol/L
end
C = round(C, 5);

model2 = model;

% Determine bounds for drug metabolites and narrow the range
model2.osenseStr = 'min';
FBA_sol = cell(length(C), 1); % Preallocate a cell array to hold the structs

% Perform Flux Balance Analysis (FBA) for each time point
tic
parfor i = 1:length(C)
    Model = model2;
    warning('off', 'all');
    Model = PBPK_constraint_WBM(Model, Physiological_param, Dosage_inf, Drug_param, C(:, i), t(i));
    warning('on', 'all');
    FBA = optimizeWBModel(Model, param);
    FBA_sol{i} = FBA;
end
toc

% Combine FBA results
FBA_sol = vertcat(FBA_sol{:});

% Save the status of the solution obtained
stat = {FBA_sol.origStat}';

% Extract time series data for fluxes
Flux_t = [];
for i = 1:length(FBA_sol)
    Flux_t(:, i) = FBA_sol(i, :).v;
end
end

