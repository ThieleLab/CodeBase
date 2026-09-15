function dX = PBPK_generic(t,X,Physiological_param,Drug_param,PBPK_model_inf,Dosage_inf)
%% Description: 
% The function produces a collection of differential equations that are
% tailored to different tissues and are parameterized for the entire body.
% Tissues: [ Adipose, Artery, Brain, Small Intestine, Large Intestine,
% Heart, Kidney, Liver, Lung, Muscle, Pancreas, Skin, Spleen, Stomach,
% Bone, Venous, Transit 1, Transit 2, Transit 3, Small intestine lumen,
% Large intestine lumen, Urine, feces]
%
% USAGE:
%   [dC] = PBPK_generic(t,C,Physiological_params,Drug_param,PBPK_model_inf,Dosage_inf)
%
% INPUT:
% t                     The simulation time point is represented by 't'. 
% C                     Concentration at a given time point
% Physiological_param  Physiological parameters Height, weight, age etc.
% Drug_param            Drug specific parameters lipophilicity, fractional
%                       unbound in plama partition coefficient.
% PBPK_model_inf        PBPK model structure and other information. 
% Dosage_inf            Dosage information administration route, ACAT
%                       parameters 
%                       'Intarvenous' 'Arterial' 'Oral' 'Intarvenous_bolus'
% dC                    Change in Concentration of the drug in a specific tissue. 
%
% OUTPUT:
% dC                Parameterized differential equations.
%
% .. Authors:
%       - Original author: Mohammad Faiz Khan.
%
% .. Last updated: 28/04/2023
%% Code 
    V = Physiological_param.V;
    Q = Physiological_param.Q;
    
    C = X(1:PBPK_model_inf.n_VTeqns);   
    C_t = X(PBPK_model_inf.n_VTeqns+1:PBPK_model_inf.n_VTeqns+PBPK_model_inf.n_EVTeqns);
    
    %  PBPK model definition: begin 
    
    Q(find(PBPK_model_inf.struct==0))=0; % Changing the flow rate value through the tissue to zero if it is not present in the study.
    
    C(find(repmat(PBPK_model_inf.struct==0,1,Drug_param.N_metabolites)))=0;
    C = C(PBPK_model_inf.G);
    dC = [];
    dC_t = [];
        for i =1:Drug_param.N_metabolites
            dc(1:23) = Organs_model([C(1+23*(i-1):23+23*(i-1));C_t],Physiological_param,Drug_param.metabolite{i},PBPK_model_inf).v;
            dc = dc.';
            dc_t = Organs_model([C(1+23*(i-1):23+23*(i-1));C_t],Physiological_param,Drug_param.metabolite{i},PBPK_model_inf).t;
            dc = diag(PBPK_model_inf.struct)*dc; % Excluding the undesired tissues that were not taken into account in the study.
            dC = [dC;dc];
            dC_t = [dC_t;dc_t];
            clear dc dc_t
        end
    %  PBPK model definition: end 
    %  PBPK model dose definition: begin 
    dA = [];
    switch Dosage_inf.Dosage_type
        case 'Arterial' % ART
                dC(2) = dC(2); %Venous blood Infusion rate
        case 'Oral' % OS
            A = X(PBPK_model_inf.n_VTeqns+PBPK_model_inf.n_EVTeqns+1:PBPK_model_inf.n_VTeqns+PBPK_model_inf.n_EVTeqns+Dosage_inf.ACAT_param.n_ACATeqns); 
            Dosage_inf.ACAT_param.C_LI = C(8); % for enterohepatic recirculation from liver define drug conc from liver.
            dA = ACAT_equ(A,Dosage_inf.ACAT_param); % dA = [dA_UND.';dA_DIS.';dA_DEG.';dA_ABS.'] 
            
            dA_ABS = dA([3*Dosage_inf.ACAT_param.n+1:4*Dosage_inf.ACAT_param.n]');
            GIR = dA_ABS(1);% Gastric abroption
            SIA = sum(dA_ABS(2:Dosage_inf.ACAT_param.n-1));% small intestine abroption
            LIA = dA_ABS(Dosage_inf.ACAT_param.n);% large intestine abroption
            
            % Connecting ACAT model and PBPK model
                dC(14) = dC(14) + GIR/Physiological_param.V(16); % gastric infusion rate 
                dC(4) = dC(4) + SIA/Physiological_param.V(4); % Small intetinal absorption  
                dC(5) = dC(5) + LIA/Physiological_param.V(5); % large intetinal absorption
        case 'Intarvenous_bolus' % IV
            input_int = (Dosage_inf.Intarvenous_bolus.dose_rate)*((t<=Dosage_inf.Intarvenous_bolus.d));
                dC(16) = dC(16) + input_int/V(16); %Venous blood Infusion rate
        case 'Intarvenous'
                dC(16) = dC(16);
        case 'EY'% EY
        case 'NA'% NA
        case 'OSIV'% OSIV
        case 'SL'% SL
        case 'Tablet'% Tablet
        case 'TD'% TD
        case 'TI'% TI
        case 'OI'% OI
        case 'OL'% OL
        case 'SC'% SC
        case 'VA'% VA
        case 'EA'% EA
        case 'IM'% IM
        case 'RE'% RE
        case 'No Dose'
    end
    %  PBPK model dose definition: end
    dC = dC(:);
    %  PBPK model metabolism definition: begin
    dM = metabolism(C(1:23*Drug_param.N_metabolites),Physiological_param,PBPK_model_inf,Drug_param.metabolite{i});
    
    dC = dC +dM;
    %  PBPK model metabolism definition: end
  % Combining all into difference vector
  if (~isempty(dC_t) && ~isempty(dA))
      dX = [dC;dC_t;dA];
  elseif (~isempty(dC_t)&& isempty(dA))
      dX = [dC;dC_t];
  elseif (isempty(dC_t)&& ~isempty(dA))
      dX = [dC;dA];
  elseif (isempty(dC_t)&& isempty(dA))
      dX = dC;
  end
end