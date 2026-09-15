function [dM] = metabolism(C,Physiological_param,PBPK_model_inf,Parameters)
% Description: 
% 
%
% USAGE:
%   [dM] = metabolism(C,Physiological_param,PBPK_model_inf,Parameters)
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
    global Sim_type
    
    if(Sim_type == 0)
        dM = PBPK_model_inf.M;
    elseif(Sim_type == 1)
        dM =eval(PBPK_model_inf.M);
    end
    
end
