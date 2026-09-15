function [Q] = organFlow(Physiological_param)
%% Description: 
% This function calulates the blood flow rate based on the physiological
% properties of the human. 
%
% USAGE:
%   [BSA] = Bodysurfacearea(Physiological_param);
%   
%
% INPUT:
% Physiological_param   Phsyiological properties such as gender height
%                       bodymass.
%
% OUTPUT:
% Q       Blood flow rate
%
% .. Authors:
%       - Original author: Mohammad Faiz Khan.
%
% .. Last updated: 28/04/2023
%% Code
    age = Physiological_param.age;
    height = Physiological_param.height;
    BSA = Physiological_param.BSA;
    C = [.01*((0.044+0.027*sex)*age+2.4*sex+3.9); % 1.adipose
         0;                                       % 2.blood
         .01*(exp(-0.48*BSA+0.04*sex+3.5));       % 3.brain   
         .01*(2*sex+14)/2;                        % 4.s int
         .01*(2*sex+14)/2;                        % 5.l int
         .01*(-0.72*height-10*sex+134);           % 6.heart
         .01*(-8.7*BSA+0.29*height-0.081*age-13); % 7.kidney
         .01*.24*(-0.108*age+1.04*sex+27.9);      % 8.liver, 24% from hep artery, rest from portal vein
         1;                                       % 9.lung
         .01*(-6.4*sex+17.5);                     %10.muscle
         .01;                                     %11.pancreas
         .05;                                     %12.skin
         .03;                                     %13.spleen
         0;                                       %14.stomach
         0;                                       %15.stomach_lumen
         0;                                       %16.sint_lumen
         0];                                      %17.lint_lumen

    C(14)=(2-sum(C)); %stomach fix

    Cardiac_output = (159*BSA-1.56*age+114);%(L/hr)
    Q = C.*Cardiac_output/60;          %15x1 vector of blood flow (L/min)

end