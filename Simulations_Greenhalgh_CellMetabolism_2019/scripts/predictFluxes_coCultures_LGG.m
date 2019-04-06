% predict fluxes for LGG with/without Caco-2 model
% only with LGG-leave

%% DMEM medium
clear all
% for adding the correct reaction names and formulas
load('rBioNetDB.mat');
initCobraToolbox
currentDir = pwd;
modelPath=currentDir;
modelIDs={
    'model_LGG'
    'model_Caco2'
    'modelJoint_Caco2_LGG'
    };

% first collect all of the reactions
Reactions={};
for i=1:length(modelIDs)
    load(strcat(modelPath,'\',modelIDs{i},'.mat'));
    Reactions=vertcat(Reactions,model.rxns);
end
Reactions=unique(Reactions);
for j=1:length(Reactions)
    MinFluxes{j+1,1}=Reactions{j,1};
    MaxFluxes{j+1,1}=Reactions{j,1};
    FluxSpans{j+1,1}=Reactions{j,1};
end
MinFluxes{1,1}='Model_Reaction_ID';
MinFluxes{1,2}='VMH_ID';
MinFluxes{1,3}='Description';
MinFluxes{1,4}='Formula';
MinFluxes{1,5}='Subsystem';
MinFluxes{1,6}='Species';
MaxFluxes{1,1}='Model_Reaction_ID';
MaxFluxes{1,2}='VMH_ID';
MaxFluxes{1,3}='Description';
MaxFluxes{1,4}='Formula';
MaxFluxes{1,5}='Subsystem';
MaxFluxes{1,6}='Species';
FluxSpans{1,1}='Model_Reaction_ID';
FluxSpans{1,2}='VMH_ID';
FluxSpans{1,3}='Description';
FluxSpans{1,4}='Formula';
FluxSpans{1,5}='Subsystem';
FluxSpans{1,6}='Species';

for i=1:length(modelIDs)
    i
    MinFluxes{1,i+6}=modelIDs{i};
    MaxFluxes{1,i+6}=modelIDs{i};
    FluxSpans{1,i+6}=modelIDs{i};
    load(strcat(modelPath,'\',modelIDs{i},'.mat'));
    % DMEM model
    model=useDMEM(model);
    % needed for Caco2 model
    model=changeRxnBounds(model,'EX_tag_hs[u]',-10,'l');
    % adjust the biomass constraints
    if strcmp(modelIDs{i},'model_LGG')
        model=changeObjective(model,'Lactobacillus_rhamnosus_GG_ATCC_53103_biomass0');
    elseif strcmp(modelIDs{i},'model_Caco2')
        model=changeObjective(model,'Host_biomass_reaction');
    end
    % set the experimentally determined growth rates as
    % constraints-microbes
    % data in Excel sheet 'bac_cellcount_beforehumix'
    if strcmp(modelIDs{i},'modelJoint_Caco2_LGG')
        model=changeRxnBounds(model,'Lactobacillus_rhamnosus_GG_ATCC_53103_biomass0',0.282377135,'b');
    end
    % now the Caco-2 cell growth rates
    % data in Excel sheet 'Bac_cellcount_Nhumix_25012017'
    if strcmp(modelIDs{i},'model_Caco2')
        model=changeRxnBounds(model,'Host_biomass_reaction',0.010632074,'b');
    end
    if strcmp(modelIDs{i},'modelJoint_Caco2_LGG')
        % no data available-here I assume 0.01 since that's the general
        % range
        model=changeRxnBounds(model,'Host_biomass_reaction',0.01,'b');
    end
    % now the SCFA measurements (only LGG alone), add as ratio
    % underlying data in Excel file 'SCFA_analysis_25_08_2017'
    % it seems there is a bug in the function so the ratio need to be
    % written in reverse for it to work.
    if strcmp(modelIDs{i},'model_LGG')
        model = addRatioReaction(model, {'EX_ac[u]' 'EX_for[u]'}, [0.418 0.326]);
    end
    % save the model
    save(strcat(modelIDs{i},'_DMEM'),'model');
    % perform computations
    formulas=printRxnFormula(model);
    currentDir = pwd;
    [minFlux,maxFlux,optsol,ret]=fastFVA(model,90,'max','ibm_cplex');
    cd(currentDir);
   for j=2:size(MinFluxes,1)
        % remove the species identifiers to get the IDs of just the reactions
        % themselves
        if strmatch('Lactobacillus_rhamnosus_GG_ATCC_53103_',MinFluxes{j,1})
            MinFluxes{j,2} = regexprep(MinFluxes{j,1},'Lactobacillus_rhamnosus_GG_ATCC_53103_','');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'\[u\]tr','\(e\)');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'IEX_','EX_');
            MinFluxes{j,6}='Lactobacillus_rhamnosus_GG_ATCC_53103';
            MaxFluxes{j,2} = MinFluxes{j,2};
            MaxFluxes{j,6}='Lactobacillus_rhamnosus_GG_ATCC_53103';
            FluxSpans{j,2} = MinFluxes{j,2};
            FluxSpans{j,6}='Lactobacillus_rhamnosus_GG_ATCC_53103';
        elseif strmatch('Host_',MinFluxes{j,1})
            MinFluxes{j,2} = regexprep(MinFluxes{j,1},'Host_','');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'\[u\]tr','\(e\)');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'\(e\)b','\(e\)');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'IEX_','EX_');
            MinFluxes{j,6}='Human';
            MaxFluxes{j,2} = MinFluxes{j,2};
            MaxFluxes{j,6}='Human';
            FluxSpans{j,2} = MinFluxes{j,2};
            FluxSpans{j,6}='Human';
        elseif strncmp('EX_',MinFluxes{j,1},3)
            MinFluxes{j,2}=MinFluxes{j,1};
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'\[u\]','\(e\)');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'glc_D','glc');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'adocbl','adpcbl');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'sbt_D','sbt-d');
            MinFluxes{j,6}='Exchange';
            MaxFluxes{j,2}=MinFluxes{j,2};
            MaxFluxes{j,6}='Exchange';
            FluxSpans{j,2}=MinFluxes{j,2};
            FluxSpans{j,6}='Exchange';
        end
        % fill in the reaction information from the database, otherwise
        % from the model
        findRxn=find(strcmp(rBioNetDB(:,1),MinFluxes{j,2}));
        if ~isempty(findRxn)
            MinFluxes{j,3}=rBioNetDB{findRxn,2};
            MinFluxes{j,4}=rBioNetDB{findRxn,3};
            MaxFluxes{j,3}=rBioNetDB{findRxn,2};
            MaxFluxes{j,4}=rBioNetDB{findRxn,3};
            FluxSpans{j,3}=rBioNetDB{findRxn,2};
            FluxSpans{j,4}=rBioNetDB{findRxn,3};
        else
            findRxn=find(strcmp(model.rxns,MinFluxes{j,1}));
            if ~isempty(findRxn)
            MinFluxes{j,3}=model.rxnNames{findRxn};
            MinFluxes{j,4}=formulas{findRxn};
            MaxFluxes{j,3}=model.rxnNames{findRxn};
            MaxFluxes{j,4}=formulas{findRxn};
            FluxSpans{j,3}=model.rxnNames{findRxn};
            FluxSpans{j,4}=formulas{findRxn};
            end
        end
        findRxn=find(strcmp(model.rxns,MinFluxes{j,1}));
        if ~isempty(findRxn)
            MinFluxes{j,5}=model.subSystems{findRxn};
            MaxFluxes{j,5}=model.subSystems{findRxn};
            FluxSpans{j,5}=model.subSystems{findRxn};
            % now fill in the fluxes
            MinFluxes{j,i+6}=minFlux(findRxn);
            MaxFluxes{j,i+6}=maxFlux(findRxn);
            if maxFlux(findRxn) > 0.0000000001 && minFlux(findRxn) > 0.0000000001
                FluxSpans{j,i+6}=maxFlux(findRxn)-minFlux(findRxn);
            elseif maxFlux(findRxn) > 0.0000000001 && minFlux(findRxn) <-0.0000000001
                FluxSpans{j,i+6}=maxFlux(findRxn) + abs(minFlux(findRxn));
            elseif maxFlux(findRxn) < -0.0000000001 && minFlux(findRxn) <-0.0000000001
                FluxSpans{j,i+6}=abs(minFlux(findRxn)) - abs(maxFlux(findRxn));
            elseif maxFlux(findRxn) > 0.0000000001 && abs(minFlux(findRxn)) <0.0000000001
                FluxSpans{j,i+6}=maxFlux(findRxn);
            elseif minFlux(findRxn) < -0.0000000001 && abs(maxFlux(findRxn)) <0.0000000001
                FluxSpans{j,i+6}=abs(minFlux(findRxn));
            elseif abs(maxFlux(findRxn)) < 0.0000000001 && abs(minFlux(findRxn)) <0.0000000001
                FluxSpans{j,i+6}=0;
            end
        else
            MinFluxes{j,i+6}=0;
            MaxFluxes{j,i+6}=0;
            FluxSpans{j,i+6}=0;
        end
    end
end
save('MinFluxes_DMEM_LGG.mat','MinFluxes');
save('MaxFluxes_DMEM_LGG.mat','MaxFluxes');
save('FluxSpans_DMEM_LGG.mat','FluxSpans');

%% SIEM medium
clear all
% for adding the correct reaction names and formulas
load('rBioNetDB.mat');
currentDir = pwd;
initCobraToolbox
cd(currentDir);
modelPath=currentDir;
modelIDs={
    'model_LGG'
    'model_Caco2'
    'modelJoint_Caco2_LGG'
    };

% first collect all of the reactions
Reactions={};
for i=1:length(modelIDs)
    load(strcat(modelPath,'\',modelIDs{i},'.mat'));
    Reactions=vertcat(Reactions,model.rxns);
end
Reactions=unique(Reactions);
for j=1:length(Reactions)
    MinFluxes{j+1,1}=Reactions{j,1};
    MaxFluxes{j+1,1}=Reactions{j,1};
    FluxSpans{j+1,1}=Reactions{j,1};
end
MinFluxes{1,1}='Model_Reaction_ID';
MinFluxes{1,2}='VMH_ID';
MinFluxes{1,3}='Description';
MinFluxes{1,4}='Formula';
MinFluxes{1,5}='Subsystem';
MinFluxes{1,6}='Species';
MaxFluxes{1,1}='Model_Reaction_ID';
MaxFluxes{1,2}='VMH_ID';
MaxFluxes{1,3}='Description';
MaxFluxes{1,4}='Formula';
MaxFluxes{1,5}='Subsystem';
MaxFluxes{1,6}='Species';
FluxSpans{1,1}='Model_Reaction_ID';
FluxSpans{1,2}='VMH_ID';
FluxSpans{1,3}='Description';
FluxSpans{1,4}='Formula';
FluxSpans{1,5}='Subsystem';
FluxSpans{1,6}='Species';

for i=1:length(modelIDs)
    MinFluxes{1,i+6}=modelIDs{i};
    MaxFluxes{1,i+6}=modelIDs{i};
    FluxSpans{1,i+6}=modelIDs{i};
    load(strcat(modelPath,'\',modelIDs{i},'.mat'));
    % SIEM medium
    model=useSIEM(model);
    % needed for Caco2/ model
    model=changeRxnBounds(model,'EX_tag_hs[u]',-10,'l');
    % adjust the biomass constraints
    if strcmp(modelIDs{i},'model_LGG')
        model=changeObjective(model,'Lactobacillus_rhamnosus_GG_ATCC_53103_biomass0');
    elseif strcmp(modelIDs{i},'model_Caco2')
        model=changeObjective(model,'Host_biomass_reaction');
    end
    % now the Caco-2 cell growth rates
    % data in Excel sheet 'Bac_cellcount_Nhumix_25012017'
    if strcmp(modelIDs{i},'model_Caco2') || strcmp(modelIDs{i},'model_Control')
        model=changeRxnBounds(model,'Host_biomass_reaction',0.00231344,'b');
    end
    if strcmp(modelIDs{i},'modelJoint_Caco2_LGG') || strcmp(modelIDs{i},'modelJoint_Control_LGG')
        model=changeRxnBounds(model,'Host_biomass_reaction',0.006414041,'b');
    end
        % save the model
    save(strcat(modelIDs{i},'_SIEM'),'model');
    % perform computations
    formulas=printRxnFormula(model);
    currentDir = pwd;
    [minFlux,maxFlux,optsol,ret]=fastFVA(model,90,'max','ibm_cplex');
    cd(currentDir);
    for j=2:size(MinFluxes,1)
        % remove the species identifiers to get the IDs of just the reactions
        % themselves
      if strmatch('Lactobacillus_rhamnosus_GG_ATCC_53103_',MinFluxes{j,1})
            MinFluxes{j,2} = regexprep(MinFluxes{j,1},'Lactobacillus_rhamnosus_GG_ATCC_53103_','');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'\[u\]tr','\(e\)');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'IEX_','EX_');
            MinFluxes{j,6}='Lactobacillus_rhamnosus_GG_ATCC_53103';
            MaxFluxes{j,2} = MinFluxes{j,2};
            MaxFluxes{j,6}='Lactobacillus_rhamnosus_GG_ATCC_53103';
            FluxSpans{j,2} = MinFluxes{j,2};
            FluxSpans{j,6}='Lactobacillus_rhamnosus_GG_ATCC_53103';
        elseif strmatch('Host_',MinFluxes{j,1})
            MinFluxes{j,2} = regexprep(MinFluxes{j,1},'Host_','');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'\[u\]tr','\(e\)');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'\(e\)b','\(e\)');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'IEX_','EX_');
            MinFluxes{j,6}='Human';
            MaxFluxes{j,2} = MinFluxes{j,2};
            MaxFluxes{j,6}='Human';
            FluxSpans{j,2} = MinFluxes{j,2};
            FluxSpans{j,6}='Human';
        elseif strncmp('EX_',MinFluxes{j,1},3)
            MinFluxes{j,2}=MinFluxes{j,1};
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'\[u\]','\(e\)');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'glc_D','glc');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'adocbl','adpcbl');
            MinFluxes{j,2}=regexprep(MinFluxes{j,2},'sbt_D','sbt-d');
            MinFluxes{j,6}='Exchange';
            MaxFluxes{j,2}=MinFluxes{j,2};
            MaxFluxes{j,6}='Exchange';
            FluxSpans{j,2}=MinFluxes{j,2};
            FluxSpans{j,6}='Exchange';
        end
        % fill in the reaction information from the database, otherwise
        % from the model
        findRxn=find(strcmp(rBioNetDB(:,1),MinFluxes{j,2}));
        if ~isempty(findRxn)
            MinFluxes{j,3}=rBioNetDB{findRxn,2};
            MinFluxes{j,4}=rBioNetDB{findRxn,3};
            MaxFluxes{j,3}=rBioNetDB{findRxn,2};
            MaxFluxes{j,4}=rBioNetDB{findRxn,3};
            FluxSpans{j,3}=rBioNetDB{findRxn,2};
            FluxSpans{j,4}=rBioNetDB{findRxn,3};
        else
            findRxn=find(strcmp(model.rxns,MinFluxes{j,1}));
            if ~isempty(findRxn)
            MinFluxes{j,3}=model.rxnNames{findRxn};
            MinFluxes{j,4}=formulas{findRxn};
            MaxFluxes{j,3}=model.rxnNames{findRxn};
            MaxFluxes{j,4}=formulas{findRxn};
            FluxSpans{j,3}=model.rxnNames{findRxn};
            FluxSpans{j,4}=formulas{findRxn};
            end
        end
        findRxn=find(strcmp(model.rxns,MinFluxes{j,1}));
        if ~isempty(findRxn)
            MinFluxes{j,5}=model.subSystems{findRxn};
            MaxFluxes{j,5}=model.subSystems{findRxn};
            FluxSpans{j,5}=model.subSystems{findRxn};
            % now fill in the fluxes
            MinFluxes{j,i+6}=minFlux(findRxn);
            MaxFluxes{j,i+6}=maxFlux(findRxn);
            if maxFlux(findRxn) > 0.0000000001 && minFlux(findRxn) > 0.0000000001
                FluxSpans{j,i+6}=maxFlux(findRxn)-minFlux(findRxn);
            elseif maxFlux(findRxn) > 0.0000000001 && minFlux(findRxn) <-0.0000000001
                FluxSpans{j,i+6}=maxFlux(findRxn) + abs(minFlux(findRxn));
            elseif maxFlux(findRxn) < -0.0000000001 && minFlux(findRxn) <-0.0000000001
                FluxSpans{j,i+6}=abs(minFlux(findRxn)) - abs(maxFlux(findRxn));
            elseif maxFlux(findRxn) > 0.0000000001 && abs(minFlux(findRxn)) <0.0000000001
                FluxSpans{j,i+6}=maxFlux(findRxn);
            elseif minFlux(findRxn) < -0.0000000001 && abs(maxFlux(findRxn)) <0.0000000001
                FluxSpans{j,i+6}=abs(minFlux(findRxn));
            elseif abs(maxFlux(findRxn)) < 0.0000000001 && abs(minFlux(findRxn)) <0.0000000001
                FluxSpans{j,i+6}=0;
            end
        else
            MinFluxes{j,i+6}=0;
            MaxFluxes{j,i+6}=0;
            FluxSpans{j,i+6}=0;
        end
    end
end
save('MinFluxes_SIEM_LGG','MinFluxes');
save('MaxFluxes_SIEM_LGG','MaxFluxes');
save('FluxSpans_SIEM_LGG','FluxSpans');


