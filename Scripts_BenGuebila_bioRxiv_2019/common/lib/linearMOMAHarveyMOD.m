function [solutionDel,totalFluxDiff,solStatus] = ...
    linearMOMAHarveyMOD(fluxDist,modelDel,osenseStr,minFluxFlag,verbFlag,rxnsList,CPLEXParamSet)
%linearMOMA Performs a linear version of the MOMA (minimization of
%metabolic adjustment) approach 
%
% [solutionDel,solutionWT,totalFluxDiff,solStatus] = 
%       linearMOMA(modelWT,modelDel,osenseStr,minFluxFlag,verbFlab)
%
%INPUTS
% modelWT           Wild type model
% modelDel          Deletion strain model
%
%OPTIONAL INPUTS
% osenseStr         Maximize ('max')/minimize ('min') (Default = 'max')
% minFluxFlag       Minimize the absolute value of fluxes in the optimal MOMA
%                   solution (Default = false)
% verbFlag          Verbose output (Default = false)
% rxnsList          indices of reacttions minimized for
% CPLEXParamSet     solver parameters file
%OUTPUTS
% solutionDel       Deletion solution structure
% solutionWT        Wild-type solution structure
% totalFluxDiff     Value of the linear MOMA objective, i.e. sum|v_wt-v_del|
% solStatus         Solution status
%
% Solves the problem
%
% min sum|v_wt - v_del|
%     S_wt*v_wt = 0
%     lb_wt <= v_wt <= ub_wt
%     c_wt'*v_wt = f_wt
%     S_del*v_del = 0
%     lb_del <= v_del <= ub_del
%
% Here f_wt is the optimal wild type objective value found by FBA
%
% Notes:
%
% 1) This formulation allows for selecting the most appropriate
% optimal wild type FBA solution as the starting point as opposed to
% picking an arbitrary starting point (original MOMA implementation).
%
% 2) The reaction sets in the two models do not have to be equal as long as
% there is at least one reaction in common
%
% Markus Herrgard 11/7/06
% Marouen Ben Guebila - added support for coupled models
%                       option to minimize for a set of rxns
                        %option to minimze to a flux

%Change Solver
changeCobraSolver('tomlab_cplex');

if (nargin <3 || isempty(osenseStr))
    osenseStr = 'max';
end
if (nargin < 4 || isempty(minFluxFlag))
    minFluxFlag = false;
end
if (nargin < 5)
    verbFlag = false;
end

% %set different reactions
% if nargin > 5
%     [indList] = setdiff(1:length(modelWT.rxns),rxnsList);
%     modelDel.rxns(indList) = cellfun(@(x) ['model1_' x],...
%         modelDel.rxns(indList),'UniformOutput', false);
% end

% LP solution tolerance
% round objective
global CBT_LP_PARAMS
if (exist('CBT_LP_PARAMS', 'var'))
    if isfield(CBT_LP_PARAMS, 'objTol')
%         tol = CBT_LP_PARAMS.objTol;
        tol = 1e-3;
    else
        tol = 1e-3;
    end
else
    tol = 1e-3;
end

fluxDist = round(fluxDist,6);

%feasibility loop
n='6';
sol.origStat=0;
while(sol.origStat~=1)
   sol=solveCobraLPCPLEX(modelDel,0,0,0,CPLEXParamSet, 0, 'ILOGcomplex');
   n=num2str(str2num(n)-1);
   CPLEXParamSet.NUMERICALEMPHASIS=1;
   CPLEXParamSet.EPMRK=0.9;
   CPLEXParamSet.ADVIND=0;
   CPLEXParamSet.EPRHS=eval(['1e-' n]);
   CPLEXParamSet.EPOPT=eval(['1e-' n]);
   if(str2num(n)==2)
       break;
   end
end


%put back the original solver parameters
CPLEXParamSet.EPRHS=1e-6;
CPLEXParamSet.EPOPT=1e-6;
% modelDel.lb(find(modelDel.c))=sol.full(find(modelDel.c));

% [nMets1,nRxns1] = size(modelWT.S);
[nMets2,nRxns2] = size(modelDel.S);
% b1 = modelWT.b;
b2 = modelDel.b;

% Match model reaction sets
% commonRxns = ismember(modelWT.rxns,modelDel.rxns);
% nCommon = sum(commonRxns);
nCommon=length(rxnsList);

if (nCommon == 0)
    error('No common rxns in the models');
end


solutionDel.f = [];
solutionDel.x = [];
solutionDel.stat = -1;


% Solve wt problem
% solutionWT = optimizeCbModel(modelWT,osenseStr);

%round solution
sol.obj = floor(sol.obj/tol)*tol;

% Variables in the following problem are
% x = [v1;v2;delta+;delta-]
% where v1 = wild type flux vector
%       v2 = deletion strain flux vector
%       delta+ = v1 - v2
%       delta- = v2 - v1


    % Construct the LHS matrix
    % Rows:
    % 1: Swt*v1 = 0 for the wild type
    % 2: Sdel*v2 = 0 for the deletion strain x
    % 3: delta+ >= v1-v2x
    % 4: delta- >= v2-v1x
    % 5: c'v1 = f1 (wild type)
%     A = [modelWT.S sparse(nMets1,nRxns2+2*nCommon);
%          sparse(nMets2,nRxns1) modelDel.S sparse(nMets2,2*nCommon);
%          createDeltaMatchMatrix(modelWT.rxns,modelDel.rxns);
%          modelWT.c' sparse(1,nRxns2+2*nCommon)];
     
    A = [modelDel.S sparse(nMets2,2*nCommon);
         createDeltaMatchMatrixMod(modelDel.rxns,rxnsList);
         modelDel.c' sparse(1,2*nCommon)];

     
    % Construct the RHS vector
    b = [b2;-fluxDist(rxnsList);fluxDist(rxnsList);sol.obj];
    
    % Construct the objective (sum of all delta+ and delta-)
    c = [zeros(nRxns2,1);ones(2*nCommon,1)];
    
    % Construct the ub/lb
    % delta+ and delta- are in [0 10000]
    lb = [modelDel.lb;zeros(2*nCommon,1)];
    ub = [modelDel.ub;10000*ones(2*nCommon,1)];

    % Construct the constraint direction vector (G for delta's, E for
    % everything else)
    if isfield(modelDel,'csense')
        csense(1:(nMets2)) = [modelDel.csense];
        csense((nMets2)+1:(nMets2+2*nCommon)+1) = 'G';
    end
    
    if (verbFlag)
        fprintf('Solving linear MOMA: %d constraints %d variables ',size(A,1),size(A,2));
    end
    
    % Solve the linearMOMA problem
    if isfield(modelDel,'csense')
        [LPproblem.A,LPproblem.b,LPproblem.c,LPproblem.lb,LPproblem.ub,...
            LPproblem.csense,LPproblem.osense] = deal(A,b,c,lb,ub,csense,1);
    else
        [LPproblem.A,LPproblem.b,LPproblem.c,LPproblem.lb,LPproblem.ub,...
            LPproblem.osense] = deal(A,b,c,lb,ub,1);
    end

    LPsolution.origStat=0;
    n='6';%feasibility tolerance
    if nargin == 7
        while(LPsolution.origStat ~= 1)%due to high load on hpc, some nmumerical 
            %instability occurs with feabile models, so we decrease the
            %feasiblity tolerance to a max of 1e-3 and stop otw.
            CPLEXParamSet
            LPsolution = solveCobraLPCPLEX(LPproblem,0,0,0,CPLEXParamSet, 0, 'ILOGcomplex');
            n=num2str(str2num(n)-1);
            CPLEXParamSet.NUMERICALEMPHASIS=1;
            CPLEXParamSet.EPMRK =0.9;
            CPLEXParamSet.ADVIND=0;
            CPLEXParamSet.EPOPT=eval(['1e-' n]);
            CPLEXParamSet.EPRHS=eval(['1e-' n]);
            if(str2num(n)==2)
                break;
            end
        end
    else
        LPsolution = solveCobraLPCPLEX(LPproblem,0);
    end

    if (verbFlag)
        fprintf('%f seconds\n',LPsolution.time);
    end
    
    if (LPsolution.origStat > 0)
        LPsolution
        solutionDel.x = LPsolution.full((1):(nRxns2));
        solutionDel.f = sum(modelDel.c.*solutionDel.x);   
        totalFluxDiff = LPsolution.obj;
    end

solutionDel.stat = LPsolution.origStat;
solStatus = LPsolution.origStat;
