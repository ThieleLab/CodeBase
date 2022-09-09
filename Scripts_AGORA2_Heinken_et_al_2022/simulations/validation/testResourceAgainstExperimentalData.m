function testResourceAgainstExperimentalData(infoFile,paths,resources,numWorkers)
% Compares the performance of one or more reconstruction resources against
% data from NJC19 (PMID:32591517), BacDive (PMID:34718743), and Madin 
% (PMID:32503990). All reconstructions in the info file with these exact
% names need to exist in all compared reconstruction resources.
%
% USAGE:
%       testResourcesAgainstExperimentalData(infoFile,paths,resources)
%
% INPUTS
%
% infoFile      Cell array with names of reconstructions to test
% paths         Cell array with paths to one or more folders with 
%               reconstructions to test
% resources     Cell array with names of reconstruction resources
%               corresponding to the folders
%
% .. Author:
% Almut Heinken, 11/2021

if length(paths) ~= length(resources)
    error('There needs to be a path to a reconstruction resource for each resource name!')
end

% save all findings
findings=struct;
for i=1:length(resources)
findings.('Uptake_NJC19').(resources{i})={};
findings.('Secretion_NJC19').(resources{i})={};
findings.('Uptake_BacDive').(resources{i})={};
findings.('Secretion_BacDive').(resources{i})={};
findings.('Enzymes_BacDive').(resources{i})={};
findings.('Uptake_Madin').(resources{i})={};
end

%% Comparison with NJC19 data

uptakeData = readInputTableForPipeline('NJC19_Uptake_Data.txt');
secretionData = readInputTableForPipeline('NJC19_Secretion_Data.txt');

[C,IA] = setdiff(uptakeData(:,1),infoFile(2:end,1),'stable');
uptakeData(IA(2:end),:)=[];
[C,IA] = setdiff(secretionData(:,1),infoFile(2:end,1),'stable');
secretionData(IA(2:end),:)=[];

TPs = zeros(size(uptakeData,1)-1,length(resources));
TNs = zeros(size(uptakeData,1)-1,length(resources));
FPs = zeros(size(uptakeData,1)-1,length(resources));
FNs = zeros(size(uptakeData,1)-1,length(resources));

tol = 0.0001;

% metabolite uptake data
for i=2:size(uptakeData,1)
    % loop through the models for the resources to compare
    for k=1:length(resources)
        model = readCbModel([paths{k} filesep uptakeData{i,1} '.mat']);
        
        % open all exchanges
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),-1000,'l');
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),1000,'u');
        
        % loop through all metabolites
        for j=2:size(uptakeData,2)
            exchange = ['EX_' uptakeData{1,j} '(e)'];
            if contains(version,'R202') % for Matlab R2020a and newer
                if uptakeData{i,j}==1
                    go=1;
                elseif uptakeData{i,j}==-1
                    go=-1;
                else
                    go=0;
                end
            else
                if strcmp(uptakeData{i,j},'1')
                    go=1;
                elseif strcmp(uptakeData{i,j},'-1')
                    go=-1;
                else
                    go=0;
                end
            end
            if go==1
                if ~isempty(find(ismember(model.rxns,exchange)))
                    model=changeObjective(model,exchange);
                    FBA=optimizeCbModel(model,'min');
                    if FBA.f < -tol
                        TPs(i-1,k) =  TPs(i-1,k) + 1;
                        findings.('Uptake_NJC19').(resources{k})(size(findings.('Uptake_NJC19').(resources{k}),1)+1,:)={uptakeData{i,1},'TP',exchange};
                    else
                        FNs(i-1,k) =  FNs(i-1,k) + 1;
                        findings.('Uptake_NJC19').(resources{k})(size(findings.('Uptake_NJC19').(resources{k}),1)+1,:)={uptakeData{i,1},'FN',exchange};
                    end
                else
                    FNs(i-1,k) =  FNs(i-1,k) + 1;
                    findings.('Uptake_NJC19').(resources{k})(size(findings.('Uptake_NJC19').(resources{k}),1)+1,:)={uptakeData{i,1},'FN',exchange};
                end
            elseif go==-1
                if ~isempty(find(ismember(model.rxns,exchange)))
                    model=changeObjective(model,exchange);
                    FBA=optimizeCbModel(model,'min');
                    if FBA.f < -tol
                        FPs(i-1,k) =  FPs(i-1,k) + 1;
                        findings.('Uptake_NJC19').(resources{k})(size(findings.('Uptake_NJC19').(resources{k}),1)+1,:)={uptakeData{i,1},'FP',exchange};
                    else
                        TNs(i-1,k) =  TNs(i-1,k) + 1;
                        findings.('Uptake_NJC19').(resources{k})(size(findings.('Uptake_NJC19').(resources{k}),1)+1,:)={uptakeData{i,1},'TN',exchange};
                    end
                else
                    TNs(i-1,k) =  TNs(i-1,k) + 1;
                    findings.('Uptake_NJC19').(resources{k})(size(findings.('Uptake_NJC19').(resources{k}),1)+1,:)={uptakeData{i,1},'TN',exchange};
                end
            end
        end
    end
end

% metabolite secretion data
for i=2:size(secretionData,1)
    % loop through the models for the resources to compare
    for k=1:length(resources)
        
        model = readCbModel([paths{k} filesep secretionData{i,1} '.mat']);
        
        % open all exchanges
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),-1000,'l');
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),1000,'u');
        
        % loop through all metabolites
        for j=2:size(secretionData,2)
            exchange = ['EX_' secretionData{1,j} '(e)'];
            if contains(version,'R202') % for Matlab R2020a and newer
                if secretionData{i,j}==1
                    go=1;
                elseif secretionData{i,j}==-1
                    go=-1;
                else
                    go=0;
                end
            else
                if strcmp(secretionData{i,j},'1')
                    go=1;
                elseif strcmp(secretionData{i,j},'-1')
                    go=-1;
                else
                    go=0;
                end
            end
            if go==1
                if ~isempty(find(ismember(model.rxns,exchange)))
                    model=changeObjective(model,exchange);
                    FBA=optimizeCbModel(model,'max');
                    if FBA.f > tol
                        TPs(i-1,k) =  TPs(i-1,k) + 1;
                        findings.('Secretion_NJC19').(resources{k})(size(findings.('Secretion_NJC19').(resources{k}),1)+1,:)={secretionData{i,1},'TP',exchange};
                    else
                        FNs(i-1,k) =  FNs(i-1,k) + 1;
                        findings.('Secretion_NJC19').(resources{k})(size(findings.('Secretion_NJC19').(resources{k}),1)+1,:)={secretionData{i,1},'FN',exchange};
                    end
                else
                    FNs(i-1,k) =  FNs(i-1,k) + 1;
                    findings.('Secretion_NJC19').(resources{k})(size(findings.('Secretion_NJC19').(resources{k}),1)+1,:)={secretionData{i,1},'FN',exchange};
                end
            elseif go==-1
                if ~isempty(find(ismember(model.rxns,exchange)))
                    model=changeObjective(model,exchange);
                    FBA=optimizeCbModel(model,'max');
                    if FBA.f > tol
                        FPs(i-1,k) =  FPs(i-1,k) + 1;
                        findings.('Secretion_NJC19').(resources{k})(size(findings.('Secretion_NJC19').(resources{k}),1)+1,:)={secretionData{i,1},'FP',exchange};
                    else
                        TNs(i-1,k) =  TNs(i-1,k) + 1;
                        findings.('Secretion_NJC19').(resources{k})(size(findings.('Secretion_NJC19').(resources{k}),1)+1,:)={secretionData{i,1},'TN',exchange};
                    end
                else
                    TNs(i-1,k) =  TNs(i-1,k) + 1;
                    findings.('Secretion_NJC19').(resources{k})(size(findings.('Secretion_NJC19').(resources{k}),1)+1,:)={secretionData{i,1},'TN',exchange};
                end
            end
        end
    end
end

% sum up the total numbers

Results=cell(1,length(resources)+1);
Results(1,2:end) = resources;

Results{2,1}='True positives';
for i=1:size(TPs,2)
    Results{2,i+1}=sum(TPs(:,i));
end
Results{3,1}='True negatives';
for i=1:size(TPs,2)
    Results{3,i+1}=sum(TNs(:,i));
end
Results{4,1}='False positives';
for i=1:size(TPs,2)
    Results{4,i+1}=sum(FPs(:,i));
end
Results{5,1}='False negatives';
for i=1:size(TPs,2)
    Results{5,i+1}=sum(FNs(:,i));
end
Results{6,1}='Sensitivity';
for i=1:size(TPs,2)
    Results{6,i+1}=(sum(TPs(:,i)))/(sum(TPs(:,i))+sum(FNs(:,i)));
end
Results{7,1}='Specificity';
for i=1:size(TPs,2)
    Results{7,i+1}=(sum(TNs(:,i)))/(sum(TNs(:,i))+sum(FPs(:,i)));
end
Results{8,1}='Accuracy';
for i=1:size(TPs,2)
    Results{8,i+1}=(sum(TPs(:,i))+sum(TNs(:,i)))/(sum(TPs(:,i))+sum(TNs(:,i))+sum(FPs(:,i))+sum(FNs(:,i)));
end
Results{9,1}='Tested models';
for i=1:size(TPs,2)
    Results{9,i+1}=size(TPs(:,i),1);
end

save([pwd filesep 'Results_NJC19'],'TPs','TNs','FPs','FNs','Results');
clear('TPs','TNs','FPs','FNs','Results');

%% Comparison with BacDive data

uptakeData = readInputTableForPipeline('BacDive_Uptake_Data.txt');
secretionData = readInputTableForPipeline('BacDive_Secretion_Data.txt');

[C,IA] = setdiff(uptakeData(:,1),infoFile(2:end,1),'stable');
uptakeData(IA(2:end),:)=[];
[C,IA] = setdiff(secretionData(:,1),infoFile(2:end,1),'stable');
secretionData(IA(2:end),:)=[];

TPs = zeros(size(uptakeData,1)-1,length(resources));
TNs = zeros(size(uptakeData,1)-1,length(resources));
FPs = zeros(size(uptakeData,1)-1,length(resources));
FNs = zeros(size(uptakeData,1)-1,length(resources));

tol = 0.0000001;

% metabolite uptake data
for i=2:size(uptakeData,1)
    % loop through the models for the resources to compare
    for k=1:length(resources)
        
        model = readCbModel([paths{k} filesep uptakeData{i,1} '.mat']);
        
        % open all exchanges
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),-1000,'l');
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),1000,'u');
        
        % loop through all metabolites
        for j=2:size(uptakeData,2)
            exchange = ['EX_' uptakeData{1,j} '(e)'];
            if contains(version,'R202') % for Matlab R2020a and newer
                if uptakeData{i,j}==1
                    go=1;
                elseif uptakeData{i,j}==-1
                    go=-1;
                else
                    go=0;
                end
            else
                if strcmp(uptakeData{i,j},'1')
                    go=1;
                elseif strcmp(uptakeData{i,j},'-1')
                    go=-1;
                else
                    go=0;
                end
            end
            if go==1
                if ~isempty(find(ismember(model.rxns,exchange)))
                    model=changeObjective(model,exchange);
                    FBA=optimizeCbModel(model,'min');
                    if FBA.f < -tol
                        TPs(i-1,k) =  TPs(i-1,k) + 1;
                        findings.('Uptake_BacDive').(resources{k})(size(findings.('Uptake_BacDive').(resources{k}),1)+1,:)={uptakeData{i,1},'TP',exchange};
                    else
                        FNs(i-1,k) =  FNs(i-1,k) + 1;
                        findings.('Uptake_BacDive').(resources{k})(size(findings.('Uptake_BacDive').(resources{k}),1)+1,:)={uptakeData{i,1},'FN',exchange};
                    end
                else
                    FNs(i-1,k) =  FNs(i-1,k) + 1;
                    findings.('Uptake_BacDive').(resources{k})(size(findings.('Uptake_BacDive').(resources{k}),1)+1,:)={uptakeData{i,1},'FN',exchange};
                end
            elseif go==-1
                if ~isempty(find(ismember(model.rxns,exchange)))
                    model=changeObjective(model,exchange);
                    FBA=optimizeCbModel(model,'min');
                    if FBA.f < -tol
                        FPs(i-1,k) =  FPs(i-1,k) + 1;
                        findings.('Uptake_BacDive').(resources{k})(size(findings.('Uptake_BacDive').(resources{k}),1)+1,:)={uptakeData{i,1},'FP',exchange};
                    else
                        TNs(i-1,k) =  TNs(i-1,k) + 1;
                        findings.('Uptake_BacDive').(resources{k})(size(findings.('Uptake_BacDive').(resources{k}),1)+1,:)={uptakeData{i,1},'TN',exchange};
                    end
                else
                    TNs(i-1,k) =  TNs(i-1,k) + 1;
                    findings.('Uptake_BacDive').(resources{k})(size(findings.('Uptake_BacDive').(resources{k}),1)+1,:)={uptakeData{i,1},'TN',exchange};
                end
            end
        end
    end
end

% metabolite secretion data
for i=2:size(secretionData,1)
    % loop through the models for the resources to compare
    for k=1:length(resources)
        
        model = readCbModel([paths{k} filesep secretionData{i,1} '.mat']);
        
        % open all exchanges
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),-1000,'l');
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),1000,'u');
        
        % loop through all metabolites
        for j=2:size(secretionData,2)
            exchange = ['EX_' secretionData{1,j} '(e)'];
            if contains(version,'R202') % for Matlab R2020a and newer
                if secretionData{i,j}==1
                    go=1;
                elseif secretionData{i,j}==-1
                    go=-1;
                else
                    go=0;
                end
            else
                if strcmp(secretionData{i,j},'1')
                    go=1;
                elseif strcmp(secretionData{i,j},'-1')
                    go=-1;
                else
                    go=0;
                end
            end
            if go==1
                if ~isempty(find(ismember(model.rxns,exchange)))
                    model=changeObjective(model,exchange);
                    FBA=optimizeCbModel(model,'max');
                    if FBA.f > tol
                        TPs(i-1,k) =  TPs(i-1,k) + 1;
                        findings.('Secretion_BacDive').(resources{k})(size(findings.('Secretion_BacDive').(resources{k}),1)+1,:)={secretionData{i,1},'TP',exchange};
                    else
                        FNs(i-1,k) =  FNs(i-1,k) + 1;
                        findings.('Secretion_BacDive').(resources{k})(size(findings.('Secretion_BacDive').(resources{k}),1)+1,:)={secretionData{i,1},'FN',exchange};
                    end
                else
                    FNs(i-1,k) =  FNs(i-1,k) + 1;
                    findings.('Secretion_BacDive').(resources{k})(size(findings.('Secretion_BacDive').(resources{k}),1)+1,:)={secretionData{i,1},'FN',exchange};
                end
            elseif go==-1
                if ~isempty(find(ismember(model.rxns,exchange)))
                    model=changeObjective(model,exchange);
                    FBA=optimizeCbModel(model,'max');
                    if FBA.f > tol
                        FPs(i-1,k) =  FPs(i-1,k) + 1;
                        findings.('Secretion_BacDive').(resources{k})(size(findings.('Secretion_BacDive').(resources{k}),1)+1,:)={secretionData{i,1},'FP',exchange};
                    else
                        TNs(i-1,k) =  TNs(i-1,k) + 1;
                        findings.('Secretion_BacDive').(resources{k})(size(findings.('Secretion_BacDive').(resources{k}),1)+1,:)={secretionData{i,1},'TN',exchange};
                    end
                else
                    TNs(i-1,k) =  TNs(i-1,k) + 1;
                    findings.('Secretion_BacDive').(resources{k})(size(findings.('Secretion_BacDive').(resources{k}),1)+1,:)={secretionData{i,1},'TN',exchange};
                end
            end
        end
    end
end

% sum up the total numbers

Results=cell(1,length(resources)+1);
Results(1,2:end) = resources;

Results{2,1}='True positives';
for i=1:size(TPs,2)
    Results{2,i+1}=sum(TPs(:,i));
end
Results{3,1}='True negatives';
for i=1:size(TPs,2)
    Results{3,i+1}=sum(TNs(:,i));
end
Results{4,1}='False positives';
for i=1:size(TPs,2)
    Results{4,i+1}=sum(FPs(:,i));
end
Results{5,1}='False negatives';
for i=1:size(TPs,2)
    Results{5,i+1}=sum(FNs(:,i));
end
Results{6,1}='Sensitivity';
for i=1:size(TPs,2)
    Results{6,i+1}=(sum(TPs(:,i)))/(sum(TPs(:,i))+sum(FNs(:,i)));
end
Results{7,1}='Specificity';
for i=1:size(TPs,2)
    Results{7,i+1}=(sum(TNs(:,i)))/(sum(TNs(:,i))+sum(FPs(:,i)));
end
Results{8,1}='Accuracy';
for i=1:size(TPs,2)
    Results{8,i+1}=(sum(TPs(:,i))+sum(TNs(:,i)))/(sum(TPs(:,i))+sum(TNs(:,i))+sum(FPs(:,i))+sum(FNs(:,i)));
end
Results{9,1}='Tested models';
for i=1:size(TPs,2)
    Results{9,i+1}=size(TPs(:,i),1);
end

save([pwd filesep 'Results_BacDive_Metabolites'],'TPs','TNs','FPs','FNs','Results');

clear('TPs','TNs','FPs','FNs','Results');

%% BacDive enzyme data

if numWorkers>0 && ~isempty(ver('parallel'))
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end

% load the data
load('BacDiveEnzymeMapping.mat')
enzymeData = readInputTableForPipeline('BacDive_Enzyme_Data.txt');

[C,IA] = setdiff(enzymeData(:,1),infoFile(2:end,1),'stable');
enzymeData(IA(2:end),:)=[];

TPs = zeros(size(enzymeData,1)-1,length(resources));
TNs = zeros(size(enzymeData,1)-1,length(resources));
FPs = zeros(size(enzymeData,1)-1,length(resources));
FNs = zeros(size(enzymeData,1)-1,length(resources));

tol = 0.0000001;

% test against enzyme data
for i=2:size(enzymeData,1)
    i
    % loop through the models for the resources to compare
    for k=1:length(resources)
        % get the right enzyme mapping
        resCol = find(strcmp(BacDiveEnzymeMapping(1,:),resources{k}));
        % read the model
        model = readCbModel([paths{k} filesep enzymeData{i,1} '.mat']);

        % open all exchanges
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),-1000,'l');
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),1000,'u');

        % loop through all enzymes
        for j=2:size(enzymeData,2)
            enzRxns = BacDiveEnzymeMapping{find(strcmp(BacDiveEnzymeMapping(:,2),enzymeData{1,j})),resCol};
            enzRxns = intersect(model.rxns,enzRxns);
            if contains(version,'R202') % for Matlab R2020a and newer
                if enzymeData{i,j}==1
                    go=1;
                elseif enzymeData{i,j}==-1
                    go=-1;
                else
                    go=0;
                end
            else
                if strcmp(enzymeData{i,j},'1')
                    go=1;
                elseif strcmp(enzymeData{i,j},'-1')
                    go=-1;
                else
                    go=0;
                end
            end
            if go==1
                if ~isempty(enzRxns)
                    % find out if any of the enzymatic reactions are active
                    % perform flux variability analysis
                    currentDir=pwd;
                    try
                        [minFlux, maxFlux, ~, ~] = fastFVA(model, 0, 'max', 'ibm_cplex', ...
                            enzRxns, 'S');
                    catch
                        warning('fastFVA could not run, so fluxVariability is instead used. Consider installing fastFVA for shorter computation times.');
                        cd(currentDir)
                        [minFlux, maxFlux] = fluxVariability(model, 0, 'max', enzRxns);
                    end
                    nonzeroRxns = enzRxns(abs(minFlux)>tol);
                    nonzeroRxns = union(nonzeroRxns,enzRxns(abs(maxFlux)>tol));
                    if ~isempty(nonzeroRxns)
                        TPs(i-1,k) =  TPs(i-1,k) + 1;
                        findings.('Enzymes_BacDive').(resources{k})(size(findings.('Enzymes_BacDive').(resources{k}),1)+1,:)={enzymeData{i,1},'TP',enzymeData{1,j}};
                    else
                        FNs(i-1,k) =  FNs(i-1,k) + 1;
                        findings.('Enzymes_BacDive').(resources{k})(size(findings.('Enzymes_BacDive').(resources{k}),1)+1,:)={enzymeData{i,1},'FN',enzymeData{1,j}};
                    end
                else
                    FNs(i-1,k) =  FNs(i-1,k) + 1;
                    findings.('Enzymes_BacDive').(resources{k})(size(findings.('Enzymes_BacDive').(resources{k}),1)+1,:)={enzymeData{i,1},'FN',enzymeData{1,j}};
                end
            elseif go==-1
                if ~isempty(enzRxns)
                    % find out if any of the enzymatic reactions are active
                    currentDir=pwd;
                    try
                        [minFlux,maxFlux] = fastFVA(model,0,'max','ibm_cplex',enzRxns);
                    catch
                        cd(currentDir)
                        [minFlux,maxFlux] = fluxVariability(model,0,'max',enzRxns);
                    end
                    nonzeroRxns = enzRxns(abs(minFlux)>tol);
                    nonzeroRxns = union(nonzeroRxns,enzRxns(abs(maxFlux)>tol));
                    if ~isempty(nonzeroRxns)
                        FPs(i-1,k) =  FPs(i-1,k) + 1;
                        findings.('Enzymes_BacDive').(resources{k})(size(findings.('Enzymes_BacDive').(resources{k}),1)+1,:)={enzymeData{i,1},'FP',enzymeData{1,j}};
                    else
                        TNs(i-1,k) =  TNs(i-1,k) + 1;
                        findings.('Enzymes_BacDive').(resources{k})(size(findings.('Enzymes_BacDive').(resources{k}),1)+1,:)={enzymeData{i,1},'TN',enzymeData{1,j}};
                    end
                else
                    TNs(i-1,k) =  TNs(i-1,k) + 1;
                    findings.('Enzymes_BacDive').(resources{k})(size(findings.('Enzymes_BacDive').(resources{k}),1)+1,:)={enzymeData{i,1},'TN',enzymeData{1,j}};
                end
            end
        end
    end
end

% sum up the total numbers

Results=cell(1,length(resources)+1);
Results(1,2:end) = resources;

Results{2,1}='True positives';
for i=1:size(TPs,2)
    Results{2,i+1}=sum(TPs(:,i));
end
Results{3,1}='True negatives';
for i=1:size(TPs,2)
    Results{3,i+1}=sum(TNs(:,i));
end
Results{4,1}='False positives';
for i=1:size(TPs,2)
    Results{4,i+1}=sum(FPs(:,i));
end
Results{5,1}='False negatives';
for i=1:size(TPs,2)
    Results{5,i+1}=sum(FNs(:,i));
end
Results{6,1}='Sensitivity';
for i=1:size(TPs,2)
    Results{6,i+1}=(sum(TPs(:,i)))/(sum(TPs(:,i))+sum(FNs(:,i)));
end
Results{7,1}='Specificity';
for i=1:size(TPs,2)
    Results{7,i+1}=(sum(TNs(:,i)))/(sum(TNs(:,i))+sum(FPs(:,i)));
end
Results{8,1}='Accuracy';
for i=1:size(TPs,2)
    Results{8,i+1}=(sum(TPs(:,i))+sum(TNs(:,i)))/(sum(TPs(:,i))+sum(TNs(:,i))+sum(FPs(:,i))+sum(FNs(:,i)));
end
Results{9,1}='Tested models';
for i=1:size(TPs,2)
    Results{9,i+1}=size(TPs(:,i),1);
end

save([pwd filesep 'Results_BacDive_Enzymes'],'TPs','TNs','FPs','FNs','Results');

clear('TPs','TNs','FPs','FNs','Results');

%% Comparison with Madin data

uptakeData = readInputTableForPipeline('Mapped_Madin_Data.txt');

[C,IA] = setdiff(uptakeData(:,1),infoFile(2:end,1),'stable');
uptakeData(IA(2:end),:)=[];

TPs = zeros(size(uptakeData,1)-1,length(resources));
FNs = zeros(size(uptakeData,1)-1,length(resources));

tol = 0.0000001;

% metabolite uptake data
for i=2:size(uptakeData,1)
    % loop through the models for the resources to compare
    for k=1:length(resources)
        
        model = readCbModel([paths{k} filesep uptakeData{i,1} '.mat']);
        
        % open all exchanges
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),-1000,'l');
        model = changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),1000,'u');
        
        % loop through all metabolites
        for j=2:size(uptakeData,2)
            exchange = ['EX_' uptakeData{1,j} '(e)'];
            if contains(version,'R202') % for Matlab R2020a and newer
                if uptakeData{i,j}==1
                    go=1;
                else
                    go=0;
                end
            else
                if strcmp(uptakeData{i,j},'1')
                    go=1;
                else
                    go=0;
                end
            end
            if go==1
                if ~isempty(find(ismember(model.rxns,exchange)))
                    model=changeObjective(model,exchange);
                    FBA=optimizeCbModel(model,'min');
                    if FBA.f < -tol
                        TPs(i-1,k) =  TPs(i-1,k) + 1;
                        findings.('Uptake_Madin').(resources{k})(size(findings.('Uptake_Madin').(resources{k}),1)+1,:)={uptakeData{i,1},'TP',exchange};
                    else
                        FNs(i-1,k) =  FNs(i-1,k) + 1;
                        findings.('Uptake_Madin').(resources{k})(size(findings.('Uptake_Madin').(resources{k}),1)+1,:)={uptakeData{i,1},'FN',exchange};
                    end
                else
                    FNs(i-1,k) =  FNs(i-1,k) + 1;
                    findings.('Uptake_Madin').(resources{k})(size(findings.('Uptake_Madin').(resources{k}),1)+1,:)={uptakeData{i,1},'FN',exchange};
                end
            end
        end
    end
end

% sum up the total numbers

Results=cell(1,length(resources)+1);
Results(1,2:end) = resources;

Results{2,1}='True positives';
for i=1:size(TPs,2)
    Results{2,i+1}=sum(TPs(:,i));
end
Results{3,1}='False negatives';
for i=1:size(TPs,2)
    Results{3,i+1}=sum(FNs(:,i));
end
Results{4,1}='Sensitivity';
for i=1:size(TPs,2)
    Results{4,i+1}=(sum(TPs(:,i)))/(sum(TPs(:,i))+sum(FNs(:,i)));
end
Results{5,1}='Tested models';
for i=1:size(TPs,2)
    Results{5,i+1}=size(TPs(:,i),1);
end

save([pwd filesep 'Results_Madin'],'TPs','FNs','Results');

save([pwd filesep 'Comparison_Findings'],'findings','-v7.3');

end
