function [dC] = Organs_model(Cint,Physiological_param,Parameters,PBPK_model_inf)
% Description: 
% This function generates set of differential equations that are parametrized for the entire body and is customized for various tissues.
% Tissues: [ Adipose, Artery, Brain, Small Intestine, Large Intestine,
% Heart, Kidney, Liver, Lung, Muscle, Pancreas, Skin, Spleen, Stomach,
% Bone, Venous, Transit 1, Transit 2, Transit 3, Small intestine lumen,
% Large intestine lumen, Urine, feces]
%
% USAGE:
%   [dC] = Organs_model(C,Physiological_param,Parameters,PBPK_model_inf)
%
% INPUT:
% dC                Concentration of the drug in a specific tissue. 
% Physiological     Q: Blood flow rates, V: Organ Volume.
% param             
% Parameters        CL_urine: Urinary clearance, ka: ,kLI ,kfecal ,CLbile ,kbile .
%
% OUTPUT:
% dC                Parameterized differential equations.
%
% .. Author:
% Mohammad Faiz Khan, Last Updated: 29/03/2023

%% Code 
    Q = Physiological_param.Q;
    V = Physiological_param.V;
    K = Parameters.K;
    CL_urine = Parameters.CL_urine;
    ka1 = Parameters.ka1;
    ka2 = Parameters.ka2;
    kLI = Parameters.kLI;
    kfecal = Parameters.kfecal;
    CL_bile = Parameters.CL_bile;
    kbile = Parameters.kbile;
    met_num = Parameters.met_num;
    % ODE_equation represents the differential equations for concentration in
    % various organs 
    % n: Number of organs 
    L=1./K; 
    % L (The purpose of writing the inverse of the Partition coefficient was to prevent the occurrence of NaN values.)
    L(isinf(L) | isnan(L)) = 0;
    % Differential equation for each tissue/compartment.
    
    C = Cint(1:length(PBPK_model_inf.struct));
    C_t = Cint(1+length(PBPK_model_inf.struct):end);
    ii = 0;
    % Adipose 
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{1+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(1),[C(1);C(2)],c_t,PBPK_model_inf.compartment_inf{1+(met_num-1)*23});
        dC.v(1) = dc.v; % vascular compartment
        dC.t{1} = dc.t;% tissue compartment
        ii = ii+length(c_t);
    % Artery
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{2+(met_num-1)*23}.n-1);
        dC.v(2) = (Q(9)*C(9)*L(9)-C(2)*Q(2))/V(2);
        ii = ii+length(c_t);
    % Brain
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{3+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(3),[C(3);C(2)],c_t,PBPK_model_inf.compartment_inf{3+(met_num-1)*23});
        dC.v(3) = dc.v;
        dC.t{3} = dc.t;
        ii = ii+length(c_t);
    % Small Intestine
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{4+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(4),[C(4);C(2)],c_t,PBPK_model_inf.compartment_inf{4+(met_num-1)*23});
        dC.v(4) = dc.v; 
        dC.t{4} = dc.t;
        ii = ii+length(c_t);
    % Large Intestine
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{5+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(5),[C(5);C(2)],c_t,PBPK_model_inf.compartment_inf{5+(met_num-1)*23});
        dC.v(5) = dc.v; 
        dC.t{5} = dc.t;
        ii = ii+length(c_t);
    % Heart
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{6+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(6),[C(6);C(2)],c_t,PBPK_model_inf.compartment_inf{6+(met_num-1)*23});
        dC.v(6) = dc.v;
        dC.t{6} = dc.t;
        ii = ii+length(c_t);
    % Kidney
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{7+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(7),[C(7);C(2)],c_t,PBPK_model_inf.compartment_inf{7+(met_num-1)*23});
        dC.v(7) = dc.v-CL_urine*C(7)*L(7)/V(7);
        dC.t{7} = dc.t; 
        ii = ii+length(c_t);
    % Liver
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{8+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(8),[C(8);C(2)],c_t,PBPK_model_inf.compartment_inf{8+(met_num-1)*23});
        dC.v(8) = dc.v+((Q(13)*C(13)*L(13)+Q(11)*C(11)*L(11)+Q(14)*C(14)*L(14)+Q(4)*C(4)*L(4)+Q(5)*C(5)*L(5))-CL_bile*C(8)*L(8)+ka1*C(20)+ka2*C(21))/V(8);
        ii = ii+length(c_t);
    % Lung
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{9+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(9),[C(9);C(16)],c_t,PBPK_model_inf.compartment_inf{9+(met_num-1)*23});
        dC.v(9) = dc.v; 
        dC.t{9} = dc.t;
        ii = ii+length(c_t);
    % Muscle
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{10+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(10),[C(10);C(2)],c_t,PBPK_model_inf.compartment_inf{10+(met_num-1)*23});
        dC.v(10) = dc.v; 
        dC.t{10} = dc.t;
        ii = ii+length(c_t);
    % Pancreas
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{11+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(11),[C(11);C(2)],c_t,PBPK_model_inf.compartment_inf{11+(met_num-1)*23});
        dC.v(11) = dc.v; 
        dC.t{11} = dc.t;
        ii = ii+length(c_t);
    % Skin
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{12+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(12),[C(12);C(2)],c_t,PBPK_model_inf.compartment_inf{12+(met_num-1)*23});
        dC.v(12) = dc.v; 
        dC.t{12} = dc.t;
        ii = ii+length(c_t);
    % Spleen
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{13+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(13),[C(13);C(2)],c_t,PBPK_model_inf.compartment_inf{13+(met_num-1)*23});
        dC.v(13) = dc.v; 
        dC.t{13} = dc.t;
        ii = ii+length(c_t);
    % Stomach
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{14+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(14),[C(14);C(2)],c_t,PBPK_model_inf.compartment_inf{14+(met_num-1)*23});
        dC.v(14) = dc.v; 
        dC.t{14} = dc.t;
        ii = ii+length(c_t);
    % Bone
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{15+(met_num-1)*23}.n-1);
        dc = compartmentEqns(Q(15),[C(15);C(2)],c_t,PBPK_model_inf.compartment_inf{15+(met_num-1)*23});
        dC.v(15) = dc.v; 
        dC.t{15} = dc.t;
        ii = ii+length(c_t);
    % Venous
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{16+(met_num-1)*23}.n-1);
        dC.v(16) = ((Q(8)*C(8)*L(8)+Q(6)*C(6)*L(6)+Q(3)*C(3)*L(3)+Q(10)*C(10)*L(10)+Q(1)*C(1)*L(1)+Q(12)*C(12)*L(12)+Q(15)*C(15)*L(15)+Q(7)*C(7)*L(7))-(Q(9)*C(16)))/V(16); 
        ii = ii+length(c_t);
    % Transit 1
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{17+(met_num-1)*23}.n-1);
        dC.v(17) = CL_bile*C(8)*L(8)-kbile*C(17); 
        ii = ii+length(c_t);
    % Transit 2
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{18+(met_num-1)*23}.n-1);
        dC.v(18) = kbile*C(17)-kbile*C(18); 
        ii = ii+length(c_t);
    % Transit 3
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{19+(met_num-1)*23}.n-1);
        dC.v(19) = kbile*C(18)-kbile*C(19); 
        ii = ii+length(c_t);
    % Small intestine lumen
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{20+(met_num-1)*23}.n-1);
        dC.v(20) = kbile*C(19)-(ka1+kLI)*C(20);
        ii = ii+length(c_t); 
    % Large intestine lumen
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{21+(met_num-1)*23}.n-1);
        dC.v(21) = kLI*C(20)-(ka2+kfecal)*C(21); 
        ii = ii+length(c_t);
    % Urine
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{22+(met_num-1)*23}.n-1);
        dC.v(22) = (CL_urine*C(7)*L(7)); 
        ii = ii+length(c_t);
    % feces
        c_t = C_t(ii+1:ii+PBPK_model_inf.compartment_inf{23+(met_num-1)*23}.n-1);
        dC.v(23) = kfecal*C(21); 
    clear L
    dC.v = dC.v(:);
    
    dC_tissue =[];
        for i=1:length(dC.t)
            dC_tissue = [dC_tissue;dC.t{i}];
        end
    clear dC.t
    dC.t = dC_tissue;
    dC.v(isnan(dC.v)|isinf(dC.v))=0;  
    dC.t(isnan(dC.t)|isinf(dC.t))=0;  
end
