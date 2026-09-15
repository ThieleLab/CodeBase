function [CO] = cardiacout(Physiological_param)
%% Description: 
% This function calulates the cardiac output based on the physiological
% properties of the human. 
%
% USAGE:
%   [CO] = cardiacout(Physiological_param);
%   
%
% INPUT:
% Physiological_param   Phsyiological properties such as gender height
%                       bodymass.
%
% OUTPUT:
% CO       cardiac output
%
%
% 
% .. Authors:
%       - Original author: Mohammad Faiz Khan.
%
% .. Last updated: 28/04/2023
%% Code
    
    age = Physiological_param.age;
    BSA = Physiological_param.BSA;
    
    if ~isfield(Physiological_param,'CO_measurement_equ')
        Physiological_param.CO_measurement_equ = 'default'; 
    end

    switch Physiological_param.CO_measurement_equ
        case 'default'
            CO = 15*(age)^0.75;
        case '1'
            CO = 159*BSA -1.56*age +114;
    end
end