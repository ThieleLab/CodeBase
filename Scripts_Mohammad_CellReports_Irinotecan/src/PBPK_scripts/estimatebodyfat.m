function [f] = estimatebodyfat(Physiological_param)
%% Description: 
% This MATLAB function calculates the percentage of body fat using the
% input variables of age, sex, height, body mass, and body surface area(BSA).
%
% USAGE:
%   [BSA] = Bodysurfacearea(Physiological_param);
%   
%
% INPUT:
% Physiological_param   Phsyiological properties such as gender height
%                       bodymass.
% OPTIONAL INPUTS:
% Ethinicity
%
% OUTPUT:
% f       percentage body fat.
%
% NOTE:
% The estimation is based on the correlation based on the paper by
%   1. Javier Gómez-Ambrosi et al (https://doi.org/10.2337%2Fdc11-1334) (for age 18 to 80)
%   2. Ernesto Cortés-Castell et al (https://doi.org/10.7717%2Fpeerj.3238) (for age 4 to 18)
% 
% .. Authors:
%       - Original author: Mohammad Faiz Khan.
%
% .. Last updated: 28/04/2023
%% Code
    BMI = Physiological_param.BMI;
    age = Physiological_param.age;
    gender = Physiological_param.gender;

    if ((4 <= age) && (age < 18))
        f = 62.627-11245.580*BMI -2 + (-259.114*BMI-1+2*310*age-0.151*age^2)*double(gender==0);
    elseif ((18<=age)&&(age<=80))
        f = -44.988 + (0.503 * age) + (10.689 * gender) + (3.172 * BMI) - (0.026 * BMI^2) + (0.181 * BMI * gender) - (0.02 * BMI * age) - (0.005 * BMI^2 * gender) + (0.00021 * BMI^2 * age);
    end
end
