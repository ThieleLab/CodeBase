function [dC] = compartmentEqns(Q,C,C_tis,compartment_inf)
%% Description: 
% -This function is designed to calculate the amount of drugs in different
%  organs, taking into account the sub-compartments within each organ.
% -The concentration calculation may vary depending on the number of assumed
%  sub-compartments within a particular organ.
% 
% USAGE:
%   [dC] = compartmentEqns(Q,C,C_tis,compartment_inf)
% 
% INPUTS:
% Q                 Blood flow rates
% V                 Organ Volume.
% L                 inverse of partition coefficient 
% C                 Concentration in the vasuclar space and plasma
% C_tis             Concentration in the extravasuclar space 
% csense            0 'Perf' 1 'Withplasma' 2 'Withoutplasma'        
% P                 Permeability  
% SA                Surface area
% n                 number of sub-compartments
% fu                fraction undbound
% HCT               Hematocrit
%
% OUTPUT:
% dC.v              Concentration of the drug in vascular compartment. 
% dC.t              Concentration of the drug in tissue compartment. 
% 
% NOTE:
% 
%
% .. Authors:
%       - Original author: Mohammad Faiz Khan.
%
% .. Last updated: May 2023
%% Code: 
    csense = compartment_inf.csense;
    HCT = compartment_inf.HCT;
    fu = compartment_inf.fu;
    n = compartment_inf.n;
    P = compartment_inf.P;
    SA = compartment_inf.SA;
    V = compartment_inf.V;
    L = compartment_inf.L;
    
    dCdt =[]; % differential equation for transport equation
    
    switch csense
        % 'Perfusion limited'
        case 0 
            dCdt = Q*(C(2)-C(1)*L(1))/V; % Basic Transport equation
            
        % 'Permeability limited With plasma'   
        case 1 
            e = zeros(n,n);
            for i=1:n-1
                A = transpose(P(i,:))*SA(i).*ones(1,2);
                A = A-2*diag(A);
                e = e + blkdiagi(A,i,n);
            end
            % Passive gradients
            E = fu*e * diag(L);
            % Final differential equation for the amount chnage in all the
            % subcompartments 
            % 'withplasma'
            disp('withplasma');
            Cgrad = [HCT*(C(2)-C(1));(1-HCT)*(C(4)-C(3));zeros(n-2,1)];
            dNdt = Q*Cgrad+E*[C(1);C(3);C_tis];
            % Assuming volume of the compartments are constant 
            dCdt = diag(V)\dNdt; % inv(diag(V))*dNdt
            
        % 'Permeability limited Without plasma compartment'    
        case 2 
            e = zeros(n,n);
            for i=1:n-1
                A = transpose(P(i,:))*SA(i).*ones(1,2);
                A = A-2*diag(A);
                e = e + blkdiagi(A,i,n);
            end
            % Passive gradients
            E = fu*e * diag(L);
            % Final differential equation for the amount change in all the
            % subcompartments 
            % 'withoutplasma'
            disp('withoutplasma');
            Cgrad = [(C(2)-C(1));zeros(n-1,1)];
            dNdt = Q*Cgrad+E*[C(1);C_tis];
            % Assuming volume of the compartments are constant 
            dCdt = diag(V)\dNdt;
    end
    
dC.v = dCdt(1); % change in concentration of drug in vascular compartment
dC.t = dCdt(2:end); % change in concentration of drug in tissue compartment
end