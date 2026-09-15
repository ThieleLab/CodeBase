function [V] = organVolM(Physiological_param)
%% Description 
% Function for determining organ volumes based on Age, Sex, Height, bodyMass, and BSA(Body Surface Area).
% It returns the V volume vector of different organs. 
% 
% INPUT:
% Physiological_param   Phsyiological properties such as gender height
%                       bodymass.
%
% OPTIONAL INPUTS:
%                * Ethinicity
% 
% OUTPUT:
% BSA       Body surface area
%
% NOTE:
%
%    The estimation is based on the correlation based on the paper by Felix
%    Stader 2019 (https://doi.org/10.1002/psp4.12399)
%    Blood weight: Maton A, Hopkins J, McLaughlin CW, Johnson S, Quon Warner M, LaHart D, Wright JD. (1993) Human Biology and Health. Englewood Cliffs NJ: Prentice Hall.
%    Kin, T., Murdoch, T. B., Shapiro, A. M. J., & Lakey, J. R. T. (2006). Estimation of Pancreas Weight from Donor Variables. Cell Transplantation, 15(2), 181–185. doi:10.3727/000000006783982133
%    
%
% .. Authors:
%       - Original author: Mohammad Faiz Khan.
%
% .. Last updated: December 2022
%% code
    age = Physiological_param.age; 
    gender = Physiological_param.gender;
    height = Physiological_param.height;
    bodyMass =  Physiological_param.bodyMass;
    BSA = Physiological_param.BSA;
    f = Physiological_param.per_bodyFat;
     m =[(0.01*f*bodyMass);                       % 1.adipose
         (0.3669 * (height*0.01)^3 + 0.03219 * bodyMass + 0.6041)*double(gender==0)+(0.3561 * (height*0.01)^3 + 0.03308 * bodyMass + 0.1833)*double(gender~=0);  % 2.blood
         exp(-0.0075*age+0.0078*height-0.9)*double(gender==0)+exp(-0.0075*age+0.0078*height-0.97)*double(gender~=0);     % 3.brain
         (0.45*3E-6*height^2.49);                    % 4.sint, 
         (0.55*3E-6*height^2.49);                    % 5.lint
         (0.34*BSA+0.0018*age-0.36);                % 6.heart
         (-0.00038*age-0.056*gender+0.33);          % 7.kidney
         exp(0.87*BSA - 0.0014*age - 1);          % 8.liver
         exp(0.028*height+0.0077*age-5.6);        % 9.lung
         (17.9*BSA-0.0667*age-5.68*gender-1.22);    %10.muscle
         ((4.355 + 0.742 * bodyMass + 0.837 * age)*double((4<=age)&&(age<18))+(-17.624 + 60.036 * BSA - 7.152 * gender)*double((18<=age)&&(age<75)))/1000;   % Pancreas                             %11.pancreas
         exp(-0.0058*age-0.37*gender+1.13);       %12.skin, dermis
         exp(1.13*BSA-3.93);                      %13.spleen
         1.05;                                    %14.stomach, not from Stader
         0;                                       %15.stomach_lumen
         0;                                       %16.small_lint_lumen
         0];                                      %17.large_lint_lumen

    %redistributing leftover weight to intestines and muscles
        fix = bodyMass - sum(m);
        m(4) = m(4)+0.1*fix;
        m(5) = m(5)+0.1*fix;
        m(10) = m(10)+0.8*fix;

     % densities (kg/L)

        p = [0.916; % 1.adipose
             1.060; % 2.blood
             1.035; % 3.brain
             1.044; % 4.gut
             1.044; % 5.gut
             1.030; % 6.heart     
             1.050; % 7.kidney     
             1.080; % 8.liver, [Heinemann et al]     
             1.050; % 9.lung     
             1.041; %10.muscle     
             1.045; %11.pancreas     
             1.116; %12.skin, dermis (epidermis & hypodermis average out to ~same)     
             1.054; %13.spleen
             1.050; %14.stomach     
             1;     %15.stomach_lumen
             1;     %16.sint_lumen
             1];    %17.lint_lumen   

    V = m./p;
end