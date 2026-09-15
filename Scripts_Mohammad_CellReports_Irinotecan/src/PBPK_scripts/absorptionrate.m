function [V] = absorptionrate(age, gender, height, bodyMass, BSA)
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
end