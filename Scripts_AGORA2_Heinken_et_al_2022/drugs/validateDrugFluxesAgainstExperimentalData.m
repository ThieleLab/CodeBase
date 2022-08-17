
%% Compute flux through drug reactions and compare against experimental data

resultsFolder = [rootDir filesep 'ComputedDrugFluxes' filesep];
mkdir(resultsFolder)

ComplexMedium = readtable('ComplexMedium.txt', 'Delimiter', '\t');
ComplexMedium=table2cell(ComplexMedium);
ComplexMedium=cellstr(string(ComplexMedium));

taxonomy = readInputTableForPipeline('AGORA2_infoFile.xlsx');

panFolder = [rootDir filesep 'panModelsAGORA2' filesep 'Species'];

inVitroData = readInputTableForPipeline('knownDrugMetabolizers.xlsx');

tol=0.0000001;

%% compute conversion potential for described cases from the literature

for i=2:size(inVitroData,1)
    % load the strain's model if applies, otherwise the pan-species model
    if isfile([refinedFolder filesep inVitroData{i,1} '.mat'])
    model=readCbModel([refinedFolder filesep inVitroData{i,1} '.mat']);
    else
        model=readCbModel([panFolder filesep 'pan' inVitroData{i,1} '.mat']);
    end
    % list exchange reactions
    exchanges = model.rxns(strncmp('EX_', model.rxns, 3));
    % open all exchanges
    model = changeRxnBounds(model, exchanges, -1000, 'l');
    model = changeRxnBounds(model, exchanges, 1000, 'u');
    % test cases with available data
    rxns=strsplit(inVitroData{i,3},',');
    rxns=intersect(rxns,model.rxns);
    if ~isempty(rxns)
        [minFlux,maxFlux]=fastFVA(model,0,'max','ibm_cplex',rxns);
        if any(minFlux>tol) || any(maxFlux>tol)
            inVitroData{i,5}='Yes';
        else
            inVitroData{i,5}='No';
        end
    else
        inVitroData{i,5}='No';
    end
end

TPs=0;
TNs=0;
FPs=0;
FNs=0;

% count results and get overview table
for i=2:size(inVitroData,1)
    if strcmp(inVitroData{i,5},'Yes')
        if strcmp(inVitroData{i,4},'Yes')
            inVitroData{i,6}='True positive';
            TPs=TPs+1;
        elseif strcmp(inVitroData{i,4},'No')
            inVitroData{i,6}='False positive';
            FPs=FPs+1;
        end
    elseif strcmp(inVitroData{i,5},'No')
        if strcmp(inVitroData{i,4},'Yes')
            inVitroData{i,6}='False negative';
            FNs=FNs+1;
        elseif strcmp(inVitroData{i,4},'No')
            inVitroData{i,6}='True negative';
            TNs=TNs+1;
        end
    end
end
writetable(cell2table(inVitroData),[resultsFolder filesep 'Overview_drug_validation.xlsx'],'FileType','spreadsheet','WriteVariableNames',false);

%%
% sensitivity or true positive rate (TPR)
TPR = TPs/(TPs + FNs);

% specificity or true negative rate (TNR)
TNR = TNs/(TNs + FPs);

% accuracy
AC = (TPs + TNs)/(TPs + TNs + FPs + FNs);
