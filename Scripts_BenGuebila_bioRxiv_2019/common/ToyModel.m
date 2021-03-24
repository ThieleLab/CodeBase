% Illustration of CRONICS approach with an Ecoli core toy mode and an ODE
% model of metabolite dynamics. This exemple computes the growth rates of
% E.coli and g6p time-course using dynamic constraints.
% Please add the COBRA Toolbox and this repository to your path
%%
% Requirements:
% Ubuntu 16.04+ / Windows 7
% MATLAB 2013b+ (https://www.mathworks.com/products/matlab.html)
% the COBRA Toolbox v3.0 (https://github.com/opencobra/cobratoolbox/)
% TOMLAB CPLEX v7.9+ or IBM CPLEX v12.6.0+
% Please check software incompatiblity at the COBRA Toolbox page (https://opencobra.github.io/cobratoolbox/latest/installation.html)
%%
initCobraToolbox
changeCobraSolver('ibm_cplex')

%%
% Case 1: Computing E.coli growth rate using indirect coupling

% Let's set the parameters
load ecoli_core_model.mat
growthVec = [];%result vector
warmStart = 0;
cpxControl.SCAIND = -1;%CPLEX parameter to control scaling infeasbilities
model.csense = repmat('E',length(model.mets),1);
% find biomass Reaction ID
growthID = findRxnIDs(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% 1. Toy ODE Model for intracellular D-glucose 6-phosphate generation over 4 hours span
% Let's set the time step (Default 1 hour to match the units of the Ecoli
% core model)
ts    = 1;
tspan = 0:ts:5;%we over extend by 1 hour to account for the last time step
y0    = 3;%initial g6p amounts in mmol/gDW
[t,y] = ode15s(@(t,y) 1*y, tspan, y0);

% We can derive two types of constraints from this model
% 1.a. generating flux of g6p (we reverse the sign to account for reaction directionality)
depg6pFlux = -1*y;

% Change model flux units based on time step (ATPM and EX_glce(e) are the only constraints in the ecoli core model)
model.lb(findRxnIDs(model,'ATPM')) = model.lb(findRxnIDs(model,'ATPM'))*ts;
model.ub(findRxnIDs(model,'ATPM')) = model.ub(findRxnIDs(model,'ATPM'))*ts;
model.lb(findRxnIDs(model,'EX_glc(e)')) = model.lb(findRxnIDs(model,'EX_glc(e)'))*ts;

% Start ODE-FBA loop
for i = 1:length(t)-1
    % Set constraints in FBA model
    % 1.b. Metabolite pooling constraints (concentration of g6p over time which is the vector y)
    metPool = y(i)-y(i+1);% change of metabolite amounts
    model.b(findMetIDs(model,'g6p[c]')) = metPool;
    % Let's assume that generating flux of g6p is the PGI glycolysis reaction
    % PGI: glucose-6-phosphate isomerase
    model = changeRxnBounds(model, 'PGI', depg6pFlux(i), 'b');
    
    % Optional step: generate solution basis for reuse in subsequent steps to
    % accelerate computation
    if warmStart==1
        if i==1
            basisReuse = 1;
        else
            basisReuse = 0;
        end
    elseif  warmStart==0
        basisReuse = 0;
    end
    
    % 2. solve FBA model for the current step step
    printLevel = 0;
    conflictResolve = 0;
    contFunctName = '';
    interface = 'ILOGcomplex';
    
    if i==1
        % For the first time step minimize the sum of fluxes (pFBA) to ensure
        % sparsity
        minNorm = 'one';
        solution = optimizeCbModel(model, 'max', minNorm);
        solution.origStat = solution.stat;
        if warmStart==1
            minNormLPCPLEX = 0;
            model.lb=solution.x;
            model.ub=solution.x;
            [solution, LPProblem] = solveCobraLPCPLEX(model, printLevel, basisReuse,...
                conflictResolve, contFunctName, minNormLPCPLEX, interface);
            solution.x=solution.full;
        end
    elseif i>1
        % For the remaining time steps:
        % Remove growth as an objective
        model.c(growthID) = 0;
        % The objective becomes minimizing the flux to the previous time
        % step (MOMA)
        try
            [solutionDel,totalFluxDiff,solStatus] = linearMOMAHarveyMOD(solution.x,model,...
                'max',1,1,1:length(model.rxns),cpxControl);
        catch ME
            solStatus    =0;
            solutionDel.x=[];
        end
        solution.origStat = solStatus; solution.x = solutionDel.x;
    end

    if solution.origStat == 1
        % Take the value of biomass from the corresponding reaction
        % and scale back the growth rate
        newGrowth=solution.x(growthID)*(1/ts);
    else
        newGrowth=0;
    end
    growthVec = [growthVec newGrowth];
end

% Plot results
hold on
plot(t(1:end-1),growthVec,'-o', 'LineWidth', 1.5, 'MarkerSize', 8);
ylabel(['growth rate (gDW/' num2str(ts) 'h)'], 'FontSize', 14)
xlabel('time (h)', 'FontSize', 14)
title('Modeling growth rate of Ecoli using indirect coupling', 'FontSize', 14)
%%
% Case 2: Computing g6p concentration time-course using direct coupling
load ecoli_core_model.mat

% Let's set the parameters
warmStart = 0;
cpxControl.SCAIND = -1;%CPLEX parameter to control scaling infeasbilities
model.csense = repmat('E',length(model.mets),1);
model.csense(findMetIDs(model,'g6p[c]')) = 'L';

% The constraints in the FBA model are not of the form 'b' like the first case
% but rather a range to allow for the correction of the ODE model.
% we add a dummy metabolite to model 0<d(g6p[c])/dt<d(g6p_ODE)/dt 
model.S(end+1,:) = model.S(findMetIDs(model,'g6p[c]'),:);
model.csense(end+1) = 'G';
model.b(end+1)=0;
model.mets{end+1}='dummy_g6p';

% ODE parameters
ts       = 0.5;%time step
timeSpan = 0:ts:4;%a ts-hour simulation setting over 5 hour period
                  % the last time step will be added in the end
y0    = 0.03;%initial g6p amounts in mmol/gDW
g6pcc = y0;
% find biomass Reaction ID
growthID = findRxnIDs(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% 1. First simulate ODE model for the current time step
% Toy ODE Model for intracellular D-glucose 6-phosphate generation over the span
% of 1 time step
tspan = [timeSpan(1) timeSpan(2)];
[t,y] = ode15s(@(t,y) 1*y, tspan, y0);
g6pccVec=[y(1) y(end)];

% Initialize constraints
% 1.a. generating flux of g6p (we reverse the sign to account for reaction directionality)
depg6pFlux = -y(1);
    
% Change model flux units based on time step (ATPM and EX_glce(e) are the only constraints in the ecoli core model)
model.lb(findRxnIDs(model,'ATPM')) = model.lb(findRxnIDs(model,'ATPM'))*ts;
model.ub(findRxnIDs(model,'ATPM')) = model.ub(findRxnIDs(model,'ATPM'))*ts;
model.lb(findRxnIDs(model,'EX_glc(e)')) = model.lb(findRxnIDs(model,'EX_glc(e)'))*ts;

% Start ODE-FBA loop
for i = 1:length(timeSpan)-1

    % Set constraints in FBA model
    % 1.b. Metabolite pooling constraints (change of metabolite amounts over time)
    metPool = g6pccVec(end) - g6pccVec(end-1);
    model.b(findMetIDs(model,'g6p[c]')) = metPool;
    % Let's assume that depleting flux of g6p is the PGI glycolysis reaction
    % PGI: glucose-6-phosphate isomerase. To illustrate the implementation
    % of intraindividual variability in the paper, we allow the reaction
    % flux to vary within a randomly generated range. In the paper, the
    % final flux values are regressed against this range to derive effect sizes
    % of each metabolic reaction.
    model = changeRxnBounds(model, 'PGI', depg6pFlux*2, 'l');
    model = changeRxnBounds(model, 'PGI', depg6pFlux*1.5, 'u');
    
    % 2.1. Optional step: generate solution basis for reuse in subsequent steps to
    % accelerate computation
    if warmStart==1
        if i==1
            basisReuse = 1;
        else
            basisReuse = 0;
        end
    elseif  warmStart==0
        basisReuse = 0;
    end
    
    % 2. Simulate FBA model
    % solve FBA model for the current step step
    printLevel = 0;
    conflictResolve = 0;
    contFunctName = '';
    interface = 'ILOGcomplex';
    
    if i==1
        % For the first time step minimize the sum of fluxes (pFBA) to ensure
        % sparsity
        minNorm = 'one';
        solution = optimizeCbModel(model, 'max', minNorm);
        solution.origStat = solution.stat;
        if warmStart==1
            minNormLPCPLEX = 0;
            model.lb=solution.x;
            model.ub=solution.x;
            [solution, LPProblem] = solveCobraLPCPLEX(model, printLevel, basisReuse,...
                conflictResolve, contFunctName, minNormLPCPLEX, interface);
            solution.x=solution.full;
        end
    elseif i>1
        % For the remaining time steps:
        % Remove growth as an objective
        model.c(growthID) = 0;
        % The objective becomes minimizing the flux to the previous time
        % step (MOMA)
        try
            [solutionDel,totalFluxDiff,solStatus] = linearMOMAHarveyMOD(solution.x,model,...
            'max',1,1,1:length(model.rxns),cpxControl);
        catch ME
            solStatus=0;
            solutionDel.x=[];
        end
        solution.origStat = solStatus; solution.x = solutionDel.x;
    end

    % 3. Solve ODE model using the new constraints from the FBA model
    if solution.origStat == 1
        depg6pFlux= solution.x(findRxnIDs(model,'PGI'))
    else
        error('The FBA problem is infeasible, please adjust the time step to 1 hour.')
    end
    % We can simply use Euler's forward method to get the new concentration
    % we also take into account the reaction directionality
    g6pccVec = [g6pccVec g6pccVec(end)-depg6pFlux];
end

% Plot results
figure;
hold on;
% we add the last time step in the simulation
plot([timeSpan timeSpan(end)+1],g6pccVec,'-bo', 'LineWidth', 1.5, 'MarkerSize', 8);
[t,y] = ode15s(@(t,y) 1*y, [0:ts:timeSpan(end)+1], y0);
plot(t,y,'-ro', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('time (h)', 'FontSize', 14)
ylabel('g6p amounts (mmol/gDW)', 'FontSize', 14)
title('Modeling g6p dynamics in E.coli using direct coupling', 'FontSize', 14);
legend('Coupled model','ODE model')