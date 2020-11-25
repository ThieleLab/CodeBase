% This script predicts known biomarker metabolites in
% different biofluid compartments (urine, blood, csf) of the whole-body
% model for inborn-errors of metabolism (IEMs).
% The meaning of the abbreviations for metabolites and IEMs used in this
% script can be found at www.vmh.life.
% The supported model options are 'male', 'female', and 'Recon3D'. Please
% define those using the variable 'sex' (e.g., sex = 'male').
%
% Ines Thiele 2018 - 2020

if ~exist('useSolveCobraLPCPLEX','var')
    global useSolveCobraLPCPLEX
    useSolveCobraLPCPLEX = 0;
end

if ~exist('resultsPath','var')
    global resultsPath
    resultsPath = which('MethodSection3.mlx');
    resultsPath = strrep(resultsPath,'MethodSection3.mlx',['Results' filesep]);
end

changeCobraSolver('ibm_cplex','LP')
changeCobraSolver('ibm_cplex','QP')
modelName = 'Harvey';

if strcmp(modelName,'Harvey')
    %load file corresponding to fileName
    male = loadPSCMfile(modelName);
    
    %standardPhysiolDefaultParameters needs to know what sex it is dealing
    %with
    sex  = male.sex;
    standardPhysiolDefaultParameters;
    
    male = physiologicalConstraintsHMDBbased(male,IndividualParameters);
    EUAverageDietNew;
    male = setDietConstraints(male, Diet);
    model = male;
    modelO = model;
end
cnt = 1;
minRxnsFluxHealthy = 1;%0.9;



%% set unified reaction constraints -- they are duplicated again in individual scripts

R = {'_ARGSL';'_GACMTRc';'_FUM';'_FUMm';'_HMR_7698';'_UAG4E';'_UDPG4E';'_GALT'; '_G6PDH2c';'_G6PDH2r';'_G6PDH2rer';...
    '_GLUTCOADHm';'_r0541'; '_ACOAD8m';'_RE2410C';'_RE2410N'};
RxnsAll2 = '';
for i = 1: length(R)
    RxnsAll = model.rxns(~cellfun(@isempty,strfind(model.rxns,R{i})));
    RxnsAll2 =[RxnsAll2;RxnsAll];
end

%excluded reactions
R2 = {'_FUMt';'_FUMAC';'_FUMS';'BBB'};
RxnsAll4 = '';
for i = 1: length(R2)
    RxnsAll3 = model.rxns(~cellfun(@isempty,strfind(model.rxns,R2{i})));
    RxnsAll4 =[RxnsAll4;RxnsAll3];
end
RxnsAll4 = unique(RxnsAll4);
IEMRxns = setdiff(RxnsAll2,RxnsAll4);
RxnMic = model.rxns(~cellfun(@isempty,strfind(model.rxns,'Micro_')));
if ~isempty(RxnMic)
    RxnMic
end
IEMRxns = setdiff(IEMRxns,RxnMic);
% set ARGSL to be irreversible
model.lb(ismember(model.rxns,IEMRxns)) = 0;

R2 = {'_r0784';'_r0463'};
RxnsAll2 = '';
for i = 1: length(R2)
    RxnsAll = model.rxns(~cellfun(@isempty,strfind(model.rxns,R2{i})));
    RxnsAll2 =[RxnsAll2;RxnsAll];
end
X = unique(RxnsAll2);
RxnMic = model.rxns(~cellfun(@isempty,strfind(model.rxns,'Micro_')));
if ~isempty(RxnMic)
    RxnMic
end
X = setdiff(X,RxnMic);
model.lb(ismember(model.rxns,X)) = 0;
model.ub(ismember(model.rxns,X)) = 0;

%%%
Rnew = {'BileDuct_EX_12dhchol[bd]_[luSI]';'BileDuct_EX_3dhcdchol[bd]_[luSI]';'BileDuct_EX_3dhchol[bd]_[luSI]';'BileDuct_EX_3dhdchol[bd]_[luSI]';'BileDuct_EX_3dhlchol[bd]_[luSI]';'BileDuct_EX_7dhcdchol[bd]_[luSI]';'BileDuct_EX_7dhchol[bd]_[luSI]';'BileDuct_EX_cdca24g[bd]_[luSI]';'BileDuct_EX_cdca3g[bd]_[luSI]';'BileDuct_EX_cholate[bd]_[luSI]';'BileDuct_EX_dca24g[bd]_[luSI]';'BileDuct_EX_dca3g[bd]_[luSI]';'BileDuct_EX_dchac[bd]_[luSI]';'BileDuct_EX_dgchol[bd]_[luSI]';'BileDuct_EX_gchola[bd]_[luSI]';'BileDuct_EX_hca24g[bd]_[luSI]';'BileDuct_EX_hca6g[bd]_[luSI]';'BileDuct_EX_hdca24g[bd]_[luSI]';'BileDuct_EX_hdca6g[bd]_[luSI]';'BileDuct_EX_hyochol[bd]_[luSI]';'BileDuct_EX_icdchol[bd]_[luSI]';'BileDuct_EX_isochol[bd]_[luSI]';'BileDuct_EX_lca24g[bd]_[luSI]';'BileDuct_EX_tchola[bd]_[luSI]';'BileDuct_EX_tdchola[bd]_[luSI]';'BileDuct_EX_tdechola[bd]_[luSI]';'BileDuct_EX_thyochol[bd]_[luSI]';'BileDuct_EX_uchol[bd]_[luSI]'};
model.ub(ismember(model.rxns,Rnew)) = 100;

modelO = model;

if 1
    model = modelO;
    
    R = {'_BUP2';'_UPPN'};
    RxnsAll = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{1}))));
    RxnsAll2 = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{2}))));
    IEMRxns = union(RxnsAll,RxnsAll2);
    RxnMic = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,'Micro_')))) ;
    IEMRxns = setdiff(IEMRxns,RxnMic);
    if ~strcmp(sex,'Recon3D')
        % add demand reactions to blood compartment for those biomarkers reported for blood
        model = addDemandReaction(model, '3uib[bc]');
        model = addDemandReaction(model, 'cala[bc]');
        model.A = model.S;
        BiomarkerRxns = {'EX_3uib[u]' 'Unknown'
            'EX_cala[u]' 'Unknown'
            'DM_3uib[bc]' 'Unknown'
            'DM_cala[bc]'  'Unknown'
            };
    else
        BiomarkerRxns = {'EX_3uib[u]' 'Unknown'
            'EX_cala[u]' 'Unknown'
            };
    end
    [IEMSol_UPB1] = checkIEM_WBM(model,IEMRxns, BiomarkerRxns,minRxnsFluxHealthy);
end



%% HAL
model = modelO;
R = '_HISD';
RxnsAll = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R))));
% exclude _HISDC reactions
RxnsAll2 = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,'_HISDC'))));
IEMRxns = setdiff(RxnsAll,RxnsAll2);
RxnMic = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,'Micro_'))));
IEMRxns = setdiff(IEMRxns,RxnMic);
if ~strcmp(sex,'Recon3D')
    % add demand reactions to blood compartment for those biomarkers reported for blood
    model = addDemandReaction(model, 'hista[bc]');
    model = addDemandReaction(model, 'his_L[bc]');
    model = addDemandReaction(model, 'urcan[bc]');
    model.A = model.S;
    BiomarkerRxns ={%'EX_hista[u]'	'Increased (blood/urine)'
        %  'DM_hista[bc]'	'Increased (blood/urine)'
        %  'EX_im4ac[u]'	'Increased (urine)'
        %  'EX_his_L[u]'	'Increased (blood/urine)'
        %  'DM_his_L[bc]'	'Increased (blood/urine)'
        'EX_urcan[u]'	'Unknown (blood/urine)'
        'DM_urcan[bc]'	'Unknown (blood/urine)'
        };
else
    BiomarkerRxns ={%'EX_hista[u]'	'Increased (blood/urine)'
        %  'EX_im4ac[u]'	'Increased (urine)'
        %  'EX_his_L[u]'	'Increased (blood/urine)'
        'EX_urcan[u]'	'Unknown (blood/urine)'
        };
end
[IEMSol_HIS] = checkIEM_WBM(model,IEMRxns, BiomarkerRxns,minRxnsFluxHealthy);

%% %% gene ID: 125061.1 % AFMID
if 0
    % nformanth is not in model
    model = modelO;
    
    R = {'_FKYNH';'_r0239';'_HMR_6717'};
    RxnsAll = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{1}))));
    RxnsAll2 = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{2}))));
    IEMRxns = union(RxnsAll,RxnsAll2);
    RxnMic = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,'Micro_')))) ;
    IEMRxns = setdiff(IEMRxns,RxnMic);
    if ~strcmp(sex,'Recon3D')
        % add demand reactions to blood compartment for those biomarkers reported for blood
        model = addDemandReaction(model, 'nformanth[bc]');
        model.A = model.S;
        BiomarkerRxns = {'EX_nformanth[u]' 'Unknown'
            'DM_nformanth[bc]' 'Unknown'
            };
    else
        BiomarkerRxns = {'EX_nformanth[u]' 'Unknown'
            };
    end
    [IEMSol_AFMID] = checkIEM_WBM(model,IEMRxns, BiomarkerRxns,minRxnsFluxHealthy);
end

%% %% gene ID: 249.1 % ALPL
if 1
    model = modelO;
    
    R = {'_r0242';'_r0587';'_r0707'};
    RxnsAll = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{1}))));
    RxnsAll2 = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{2}))));
    IEMRxns = union(RxnsAll,RxnsAll2);
    RxnMic = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,'Micro_')))) ;
    IEMRxns = setdiff(IEMRxns,RxnMic);
    if ~strcmp(sex,'Recon3D')
        % add demand reactions to blood compartment for those biomarkers reported for blood
        model = addDemandReaction(model, 'ethamp[bc]');
        model.A = model.S;
        BiomarkerRxns = {'EX_ethamp[u]' 'Unknown'
            'DM_ethamp[bc]' 'Unknown'
            };
    else
        BiomarkerRxns = {'EX_ethamp[u]' 'Unknown'
            };
    end
    [IEMSol_ALPL] = checkIEM_WBM(model,IEMRxns, BiomarkerRxns,minRxnsFluxHealthy);
end

%% %% gene ID: 29958.1 % DMGDH
if 1
    model = modelO;
    
    R = {'_DMGDHm'};
    RxnsAll2 = '';
    for i = 1: length(R)
        RxnsAll = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{i}))));
        RxnsAll2 =[RxnsAll2;RxnsAll];
    end
    IEMRxns = unique(RxnsAll2);
    RxnMic = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,'Micro_')))) ;
    IEMRxns = setdiff(IEMRxns,RxnMic);
    if ~strcmp(sex,'Recon3D')
        % add demand reactions to blood compartment for those biomarkers reported for blood
        model = addDemandReaction(model, 'dmgly[bc]');
        model.A = model.S;
        BiomarkerRxns = {'EX_dmgly[u]' 'Unknown'
            'DM_dmgly[bc]' 'Unknown'
            };
    else
        BiomarkerRxns = {'EX_dmgly[u]' 'Unknown'
            };
    end
    [IEMSol_DMGDH] = checkIEM_WBM(model,IEMRxns, BiomarkerRxns,minRxnsFluxHealthy);
end

if 0
    model = modelO;
    
    R = {'_DMGDHm';'_HMR_4700';'_HMR_4701'};
    RxnsAll = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{1}))));
    RxnsAll2 = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{2}))));
    IEMRxns = union(RxnsAll,RxnsAll2);
    RxnMic = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,'Micro_')))) ;
    IEMRxns = setdiff(IEMRxns,RxnMic);
    if ~strcmp(sex,'Recon3D')
        model = addDemandReaction(model, 'dmgly[bc]');
        model.A = model.S;
        BiomarkerRxns = {'EX_dmgly[u]' 'Unknown'
            'DM_dmgly[bc]' 'Unknown'
            };
    else
        BiomarkerRxns = {'EX_dmgly[u]' 'Unknown'
            };
    end
    [IEMSol_DMGDH2] = checkIEM_WBM(model,IEMRxns, BiomarkerRxns,minRxnsFluxHealthy);
end


if 1
    model = modelO;
    
    R = {'_MMSAD1m';'_MMSAD3m'};
    RxnsAll = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{1}))));
    RxnsAll2 = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{2}))));
    IEMRxns = union(RxnsAll,RxnsAll2);
    RxnMic = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,'Micro_')))) ;
    IEMRxns = setdiff(IEMRxns,RxnMic);
    if ~strcmp(sex,'Recon3D')
        model = addDemandReaction(model, '3hmp[bc]');
        model.A = model.S;
        BiomarkerRxns = {'EX_3hmp[u]' 'Unknown'
            'DM_3hmp[bc]' 'Unknown'
            };
    else
        BiomarkerRxns = {'EX_3hmp[u]' 'Unknown'
            };
    end
    [IEMSol_ALDH6A1] = checkIEM_WBM(model,IEMRxns, BiomarkerRxns,minRxnsFluxHealthy);
end

if 1
    model = modelO;
    
    R = {'_ACODA';'_RE2640C'};
    RxnsAll = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{1}))));
    RxnsAll2 = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,R{2}))));
    IEMRxns = union(RxnsAll,RxnsAll2);
    RxnMic = model.rxns(find(~cellfun(@isempty,strfind(model.rxns,'Micro_')))) ;
    IEMRxns = setdiff(IEMRxns,RxnMic);
    if ~strcmp(sex,'Recon3D')
        model = addDemandReaction(model, 'acgly[bc]');
        model = addDemandReaction(model, 'acglu[bc]');
        model.A = model.S;
        BiomarkerRxns = {'EX_acgly[u]' 'Unknown'
            'DM_acgly[bc]' 'Unknown'
            'EX_acglu[u]' 'Unknown'
            'DM_acglu[bc]' 'Unknown'
            };
    else
        BiomarkerRxns = {'EX_acgly[u]' 'Unknown'
            'EX_acglu[u]' 'Unknown'
            };
    end
    [IEMSol_ACY1] = checkIEM_WBM(model,IEMRxns, BiomarkerRxns,minRxnsFluxHealthy);
end

clear Bio* Do* H_* IEM IEMRxns R R2 RxnsA* Un* Up* X cnt i j minR* model vars*
%clearvars -except Table_sum_results Accuracy Precision FalseDiscoveryRate NumDiseases NumBiomarkers


vars = who;
vars_IEM = strmatch('IEMSol_',vars);

clear Table_IEM
cnt =1;
for i = 1 : length(vars_IEM)
    % read in IEM solutions
    clear IEM
    IEM = evalin('base',vars{vars_IEM(i)});
    
    for j = 5 : 2 : size(IEM,1)
        
        % create new table with all results
        Table_IEM{cnt,1} = regexprep(vars{vars_IEM(i)},'IEMSol_',''); % IEM abbr
        Table_IEM{cnt,2} = regexprep(IEM{j,1},'Healthy:',''); % biomaker
        Table_IEM{cnt,3} = (IEM(j,2)); % healthy original values
        Table_IEM{cnt,4} = (IEM(j+1,2)); % disease original values
        
        cnt = cnt +1;
    end
end

if strcmp(modelName,'Harvey')
  
        save(['Results_IEM_Harvey_1_03_SimulationResults'])

end
