%printLevel
print = 1;
%relax FBA options
relaxOption.internalRelax     = 2;
relaxOption.exchangeRelax     = 0;
relaxOption.steadyStateRelax  = 0;
%Optimize model
harvey.A = harvey.S;
changeCobraSolver('tomlab_cplex','LP');
harvey.csense = [harvey.csense;repmat('E',length(harvey.mets)-length(harvey.csense),1)];
if isunix
    cpxControl = CPLEXParamSet_Harvey;
else
    cpxControl.SCAIND = -1;
end
if isequal(simType,'pFBA')
    FBA=solveCobraLPCPLEX(harvey,print,0,0,cpxControl);
    if ~FBA.stat
        %relax infeasible models
        fprintf('sol not found \n')
        sol = relaxedFBA(harvey,relaxOption,'l1',cpxControl);
        relaxedRxn=unique([find(sol.p);find(sol.q)]);
        fprintf('There are %d relaxed constraints \n',...
                length(relaxedRxn));
         harvey.lb=harvey.lb-sol.p;
         harvey.ub=harvey.ub+sol.q;
    end
    FBA=minFluxSum(harvey,0,cpxControl,print);
elseif isequal(simType,'FBA')
    FBA=solveCobraLPCPLEX(harvey,print,0,0,cpxControl);
    if ~FBA.stat
        %relax infeasible models
        fprintf('sol not found \n')
        sol = relaxedFBA(harvey,relaxOption,'l1',cpxControl);
        relaxedRxn=unique([find(sol.p);find(sol.q)]);
        fprintf('There are %d relaxed constraints \n',...
                length(relaxedRxn));
        harvey.lb=harvey.lb-sol.p;
        harvey.ub=harvey.ub+sol.q;
        FBA=solveCobraLPCPLEX(harvey,print,0,0,cpxControl,0);
    end
elseif isequal(simType,'CRONICS')
    if i==1
        %check model feasibility
        FBA=solveCobraLPCPLEX(harvey,print,0,0,cpxControl);
        if ~FBA.stat
            %relax infeasible models
            fprintf('sol not found \n')
            sol = relaxedFBA(harvey,relaxOption,'l1',cpxControl);
            relaxedRxn=unique([find(sol.p);find(sol.q)]);
            fprintf('There are %d relaxed constraints \n',...
                length(relaxedRxn));
            harvey.lb=harvey.lb-sol.p;
            harvey.ub=harvey.ub+sol.q;
        end
        FBA=minFluxSum(harvey,0,cpxControl,print);
    else
        FBA=solveCobraLPCPLEX(harvey,print,0,0,cpxControl);
        if ~FBA.stat
            %relax infeasible models
            fprintf('sol not found \n')
            sol = relaxedFBA(harvey,relaxOption,'l1',cpxControl);
            relaxedRxn=unique([find(sol.p);find(sol.q)]);
            fprintf('There are %d relaxed constraints \n',...
               length(relaxedRxn));

            harvey.lb=harvey.lb-sol.p;
            harvey.ub=harvey.ub+sol.q;
        end
        rxnsList = [...
            find(cellfun(@(x) isequal(x,'Kidney'),harvey.organs))'...
            find(cellfun(@(x) isequal(x,'Pancreas'),harvey.organs))'...
            find(cellfun(@(x) isequal(x,'Brain'),harvey.organs))'...
            find(cellfun(@(x) isequal(x,'Muscle'),harvey.organs))'...
            find(cellfun(@(x) isequal(x,'Adipocytes'),harvey.organs))'...
            find(cellfun(@(x) isequal(x,'Liver'),harvey.organs))'...
            find(cellfun(@(x) isequal(x,'NA'),harvey.organs))' ...
            find(harvey.c)'];%plasma liver and obj fun
        [solutionStruc,totalFluxDiff,solStatus] = ...
          linearMOMAHarveyMOD(prevTSFD,harvey,'max',0,1,rxnsList,cpxControl);
        FBA.stat=solStatus;
        FBA.obj =solutionStruc.x(find(harvey.c));
        FBA.full=solutionStruc.x;
        sprintf('status is %d and obj is %f',FBA.stat,FBA.obj)
    end
end