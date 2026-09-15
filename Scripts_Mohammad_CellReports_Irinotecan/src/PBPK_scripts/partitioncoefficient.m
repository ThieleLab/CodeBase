function [K] = partitioncoefficient(Params,Method)
%% Description: 
% This is a MATLAB function that calculates the partition coefficient of
% drugs in pharmacokinetics. It generates a vector of K values representing
% the ratio of drug concentration in tissue to that in blood for various
% organs. Please note that if you require the partition coefficient between
% tissue and plasma (K_(T:P)), you can use the relation K_(T:B)=K_(T:P)/BP. 
% 
% USAGE:
%   [K] = partitioncoefficient(Params,Method).
%   
% INPUTS:
%     Params
%         logP              lipophilicity.
%         f_up              fractional unbound plasma.
%         BP                Blood to plasma rartio
%         pKa               dissociation coefficient in case of Diprotic Acid give
%                           a vector of [pKa1,pKa2], In these equations pKa1<pKa2. 
%         pKb               dissociation coefficient in case of Diprotic base give
%                           a vector of [pKb1,pKb2], In these equations pKb1<pKb2 
%         logD              if logD is not given use logP.
%     Method            'Poulin', 'Rodger', 'Arundel'
% 
% OPTIONAL INPUTS:
%     varargin          nature of compounds in case of rodger's method of
%                       partition coefficient prediction.
%                       'Acid','Base','Diprotic_acid','Diprotic_base',
%                       'Weak_Diprotic_base','Neutral','Zwitterions_sp',
%                       'Zwitterions'.
%                   
% Output: 
%     K                 partition coefficient vector (Tissue to blood)
% 
% NOTE:
%  The estimation is based on the correlation based on the paper by
%  1. Physiologically Based Pharmacokinetic Modeling 1: Predicting the
%     Tissue Distribution of Moderate-to-Strong Bases	https://doi.org/10.1002/jps.20322.
%  2. Prediction of Pharmacokinetics Prior to In Vivo Studies. 1.
%     Mechanism?Based Prediction of Volume of Distribution	https://doi.org/10.1002/jps.10005.
%  3. Prediction of Prediction of pharmacokinetics prior studies. II.
%     Generic physiologically based pharmacokinetic models of drug disposition	https://doi.org/10.1002/jps.10128.
%  4. Rodgers T, Rowland M. 2006. Physiologically?based Pharmacokinetic
%     Modeling 2: Predicting the tissue distribution of acids, very weak
%     bases, neutrals and zwitterions. J Pharm Sci 95:1238–1257.	https://doi.org/10.1002/jps.20857.
%  5. Ye, M., Nagar, S., and Korzekwa, K. (2016) A physiologically based
%     pharmacokinetic model to predict the pharmacokinetics of highly
%     protein-bound drugs and the impact of errors in plasma protein binding.
%     Biopharm. Drug Dispos., 37: 123– 141. doi: 10.1002/bdd.1996.   
%
% .. Authors:
%       - Original author: Mohammad Faiz Khan.
%
% .. Last updated: 28/04/2023
%% Code: 
logP = Params.logP;
logD = Params.logD;
f_up = Params.f_up;
BP = Params.BP;
pKa = Params.pKa;
pKb = Params.pKb;

switch Method
    case 'Poulin'
        %% Poulin's Method for species specific partition coefficient prediction.
        %   for Human 
        %    V_nlp  V_php  V_water
        CK =[.79    .002   .180; % 1.Adipose
            .0035  .00225 .945; % 2.Artery (blood/plasma)
            .051   .0565  .770; % 3.brain
            .0487  .0163  .718; % 4.SInt
            .0487  .0163  .718; % 5.LInt
            .0115  .0166  .758; % 6.heart
            .0207  .0162  .783; % 7.kidney     
            .0348  .0252  .751; % 8.liver     
            .003   .009   .811; % 9.lung     
            .0238  .0072  .760; %10.muscle
            .0723  .0188  .660; %11.pancreas 
            .0284  .0111  .718; %12.skin     
            .0201  .0198  .788; %13.spleen      
            .0338  .0182  .784; %14.stomach 
            0.074 0.0011 0.439; %15.Bone 
            .0035  .00225 .945; % 16.Venous(blood/plasma)
            0       0       0]; %17.lint_lumen
        % for Rat     
        % CK=[ 0.853    0.002   0.12 ;
        %      0.00147  0.00083 0.96;
        %      0.0392   0.0533  0.788 ;
        %      0.0292   0.0138  0.749;
        %      0.0292   0.0138  0.749;
        %      0.014    0.0118  0.779;
        %      0.0123   0.0284  0.771;
        %      0.0138   0.0303  0.705 ;
        %      0.0219   0.014   0.79 ;
        %      0.01     0.009   0.756 ;
        %      .0723    .0188   .660;
        %      0.0239   0.018   0.651 ;
        %      0.0077   0.0136  0.771 ;
        %      .0338    .0182   .784;
        %      0.0273   0.0027  0.446 ;
        %      0.00147  0.00083 0.96;
        %      0        0       0]; 
        
        POW = 10^logP;
        DOW = 10^logD;
        fut=1/(1+((1-f_up)*0.5/(f_up)));
        % non adipose
        K(:,1)=((POW*(CK(:,1)+.3*CK(:,2))+CK(:,3)+.7*CK(:,2))/(POW*(CK(2,1)+0.3*CK(2,2))+CK(2,3)+.7*CK(2,2)))*(f_up/fut);
        % adipose, assume logD ~ logP and fut = 1 for adipose due to minimal binding proteins
        K(1,1)=((DOW*(CK(1,1)+.3*CK(1,2))+CK(1,3)+.7*CK(1,2))/(DOW*(CK(2,1)+0.3*CK(2,2))+CK(2,3)+.7*CK(2,2)))*(f_up/1);
        K = K/BP;   % Converting from T:P(tissue to plasma) to T:B
    case 'Rodger'
        natures_cmpd = Params.natures_cmpd;
    % Rodger and Rowland Method for species specific partition coefficient prediction.    
        pHiw=7.4; %https://doi.org/10.5599%2Fadmet.638
        pHp=7; %https://doi.org/10.5599%2Fadmet.638
        pHbc=7.22; %https://doi.org/10.5599%2Fadmet.638
        Haematocrit= 0.44; %https://doi.org/10.5599%2Fadmet.638
        P1 =10^(logP);
        P2 =10^(logD);
        fNLp = 0.0023;
        fNPp = 0.0013;
        % fNPlipids flipids fEW     fIW     fAR     fLR     CAP (Taken from https://doi.org/10.1002/jps.20322  and 10.1002/bdd.1996)
        % R= [nph       nl      ew      iw      ar      lr      ap];
        % Same for rat and human 
        R=[0.0016	0.853   0.017	0.135	0.049	0.068   0.4;    % 1.Adipose
           0.0029	0.0017	0       0.603	0       0       0.5;    % 2.Artery (blood/plasma)
           0.0015	0.039	0.162	0.62	0.048	0.041   0.4;    % 3.brain
           0.0125	0.038	0.282	0.475	0.158	0.141   2.41;   % 4.SInt
           0.0125	0.038	0.282	0.475	0.158	0.141   2.41;   % 5.LInt
           0.0111	0.014	0.32	0.475	0.157	0.16    2.25;   % 6.heart
           0.0242	0.012	0.273	0.483	0.13	0.137   5.03;   % 7.kidney     
           0.024	0.014	0.161	0.573	0.086	0.161   4.56;   % 8.liver     
           0.0128	0.022	0.336	0.446	0.212	0.168   3.91;   % 9.lung     
           0.0072	0.01	0.118	0.63	0.064	0.059   1.53;   %10.muscle
           0.0093	0.0403	0.12	0.664	0.06	0.06    1.67;   %11.pancreas 
           0.0044	0.06	0.382	0.291	0.277	0.096   1.32;   %12.skin     
           0.0113	0.0077	0.207	0.579	0.097	0.207   3.18;   %13.spleen      
           0        0       0       0       0       0       0;      %14.stomach 
           0.0017	0.017	0.1     0.346	0.1     0.05    0.67;   %15.Bone 
           0.0029	0.0017	0       0.603	0       0       0.5;    % 16.Venous(blood/plasma)
           0        0       0       0       0       0       0];     %17.lint_lumen
           switch natures_cmpd
            case 'Acid'
                % For acids 
                X=1+10^(pHiw-pKa(1));
                Y=1+10^(pHp-pKa(1));
                % For rest of the tissue
                Kpu(:,1)=R(:,3)+X*R(:,4)/Y+((P1*R(:,2)+(0.3*P1+0.7)*R(:,1))/Y)+(1/f_up-1-(P1*fNLp+(0.3*P1+0.7)*fNPp)/Y)*R(:,5);
                % For adipose tissue
                Kpu(1,1)=R(1,3)+X*R(1,4)/Y+((P2*R(1,2)+(0.3*P2+0.7)*R(1,1))/Y)+(1/f_up-1-(P2*fNLp+(0.3*P2+0.7)*fNPp)/Y)*R(1,5);     
                K = Kpu*f_up/BP; 
            case 'Base'
                % For base
                X=1+10^(pKb(1)-pHiw);    
                Y=1+10^(pKb(1)-pHp);
                X1=1+10^(pKb(1)-pHbc);    
                Y1=1+10^(pKb(1)-pHp);
                X2=10^(pKb(1)-pHbc);    
                KpuBC=(BP-1+Haematocrit)/Haematocrit/f_up;
                % For rest of the tissue
                KaAP=(KpuBC-(X1/Y1)*R(2,4)-(P1*R(2,2)+(0.3*P1+0.7)*R(2,1))/Y1)*(Y1/R(2,7)/X2);
                Kpu(:,1)=R(:,3)+X*R(:,4)/Y+(P1*R(:,2)+(0.3*P1+0.7)*R(:,1))/Y+(KaAP*R(:,7)*(X-1))/Y;
                % For adipose tissue
                KaAP=(KpuBC-(X1/Y1)*R(2,4)-(P2*R(2,2)+(0.3*P2+0.7)*R(2,1))/Y1)*(Y1/R(2,7)/X2);
                Kpu(1,1)=R(1,3)+X*R(1,4)/Y+(P2*R(1,2)+(0.3*P2+0.7)*R(1,1))/Y+(KaAP*R(1,7)*(X-1))/Y;
                K = Kpu*f_up/BP; 
            case 'Diprotic_acid'
                % Diprotic acids
                % In this equations pKa1<pKa2 
                X=1+10^(pHiw-pKa(1))+10^(-pKa(2)-pKa(1)+2*pHiw);
                Y=1+10^(pHp-pKa(1))+10^(-pKa(2)-pKa(1)+2*pHp);
                % For rest of the tissue
                Kpu(:,1)= R(:,3) + X*R(:,4)/Y + ((P1*R(:,2)+(0.3*P1+0.7)*R(:,1))/Y) + (1/f_up-1-(P1*fNLp+(0.3*P1+0.7)*fNPp)/Y)*R(:,7); 
                % For adipose tissue
                Kpu(1,1)= R(1,3) + X*R(1,4)/Y + ((P1*R(1,2)+(0.3*P1+0.7)*R(1,1))/Y) + (1/f_up-1-(P1*fNLp+(0.3*P1+0.7)*fNPp)/Y)*R(1,7); 
                K = Kpu*f_up/BP; 
            case 'Diprotic_base'
                % Diprotic bases
                % In these equations pKb1<pKb2 
                X=1+10^(pKb(2)-pHiw)+10^(pKb(2)+pKb(1)-2*pHiw);
                Y=1+10^(pKb(2)-pHp)+10^(pKb(2)+pKb(1)-2*pHp);
                X1=1+10^(pKb(2)-pHbc)+10^(pKb(2)+pKb(1)-2*pHbc);
                X2=10^(pKb(2)-pHbc)+10^(pKb(2)+pKb(1)-2*pHbc);
                Y1=1+10^(pKb(2)-pHp)+10^(pKb(2)+pKb(1)-2*pHp);
                
                KpuBC = (BP-1+haematocrit)/haematocrit/f_up;
                % For rest of the tissue
                KaAP = (KpuBC-(X1/Y1)*R(2,4)-((P1*R(2,2)+(0.3*P1+0.7)*R(2,1))/Y1))*(Y1/R(2,7)/X2);
                Kpu(:,1) = R(:,3)+X*R(:,4)/Y+(P1*R(:,2)+(0.3*P1+0.7)*R(:,1))/Y+(KaAP*R(:,7)*(X-1))/Y;
                % For adipose tissue
                KaAP = (KpuBC-(X1/Y1)*R(2,4)-((P2*R(2,2)+(0.3*P2+0.7)*R(2,1))/Y1))*(Y1/R(2,7)/X2);
                Kpu(1,1) = R(:,3)+X*R(:,4)/Y+(P2*R(:,2)+(0.3*P2+0.7)*R(:,1))/Y+(KaAP*R(:,7)*(X-1))/Y;
                K = Kpu*f_up/BP; 
            case 'Weak_Diprotic_base'
                % For very weak Diprotic base
                % In these equations pKb1<pKb2 
                X = 1+10^(pKb(2)-pHiw)+10^(pKb(2)+pKb(1)-2*pHiw);
                Y = 1+10^(pKb(2)-pHp)+10^(pKb(2)+pKb(1)-2*pHp);
                % R= [nph       nl      ew      iw      ar      lr      ap];
                % For rest of the tissue
                Kpu(:,1) = R(:,3)+X*R(:,4)/Y+((P1*R(:,2)+(0.3*P1+0.7)*R(:,1))/Y)+(1/f_up-1-(P1*fNLp+(0.3*P1+0.7)*fNPp)/Y)*R(:,5);
                % For adipose tissue
                Kpu(1,1) = R(1,3)+X*R(1,4)/Y+((P2*R(1,2)+(0.3*P2+0.7)*R(1,1))/Y)+(1/f_up-1-(P2*fNLp+(0.3*P2+0.7)*fNPp)/Y)*R(1,5);
                K = Kpu*f_up/BP; 
            case 'Neutral'
                % Neutrals
                X=1;
                Y=1; 
                % For rest of the tissue
                Kpu(:,1) = X*R(:,4)/Y+R(:,3)+((P1*R(:,2)+(0.3*P1+0.7)*R(:,1))/Y)+(1/f_up-1-(P1*fNLp+(0.3*P1+0.7)*fNPp)/Y)*R(:,6);
                % For adipose tissue
                Kpu(1,1) = X*R(1,4)/Y+R(1,3)+((P2*R(1,2)+(0.3*P2+0.7)*R(1,1))/Y)+(1/f_up-1-(P2*fNLp+(0.3*P2+0.7)*fNPp)/Y)*R(1,6);
                K = Kpu*f_up/BP; 
            case 'Zwitterions_sp'
                % For Zwitterions, with at least one basic pKa>7
                X=1+10^(pKb-pHiw)+10^(pHiw-pKa);
                Y=1+10^(pKb-pHp)+10^(pHp-pKa);
                X1=1+10^(pKb-pHbc)+10^(pHbc-pKa); 
                Y1=1+10^(pKb-pHp)+10^(pHp-pKa);
                X2=10^(pKb-pHbc)+10^(pHbc-pKa);
                KpuBC=(BP-1+haematocrit)/haematocrit/f_up; 
                % R= [nph       nl      ew      iw      ar      lr      ap];
                % For rest of the tissue
                KaAP=(KpuBC-(X1/Y1)*R(2,4)-(P1*R(2,2)+(0.3*P1+0.7)*R(2,1))/Y1)*(Y1/R(2,7)/X2); 
                Kpu(:,1)=R(:,3)+X*R(:,4)/Y+(P1*R(:,2)+(0.3*P1+0.7)*R(:,1))/Y+((KaAP*R(:,7)*10^(pKb-pHiw))+10^(pHiw-pKa))/Y;
                % For adipose tissue
                KaAP=(KpuBC-(X1/Y1)*R(2,4)-(P1*R(2,2)+(0.3*P1+0.7)*R(2,1))/Y1)*(Y1/R(2,7)/X2); 
                Kpu(1,1)=R(1,3)+X*R(1,4)/Y+(P2*R(1,2)+(0.3*P2+0.7)*R(1,1))/Y+((KaAP*R(1,7)*10^(pKb-pHiw))+10^(pHiw-pKa))/Y;
                
                K = Kpu*f_up/BP; 
            case 'Zwitterions'
                % For All other zwitterions
                X = 1+10^(pKb-pHiw)+10^(pHiw-pKa);
                Y = 1+10^(pKb-pHp)+10^(pHp-pKa);
                % R= [nph       nl      ew      iw      ar      lr      ap];
                % For rest of the tissue
                Kpu(:,1) = R(:,3) + X*R(:,4)/Y+((P1*R(:,2)+(0.3*P1+0.7)*R(:,1))/Y)+(1/f_up-1-(P1*fNLp+(0.3*P1+0.7)*fNPp)/Y)*R(:,5); 
                % For adipose tissue
                Kpu(1,1) = R(1,3) + X*R(1,4)/Y+((P2*R(1,2)+(0.3*P2+0.7)*R(1,1))/Y)+(1/f_up-1-(P2*fNLp+(0.3*P2+0.7)*fNPp)/Y)*R(1,5);
                K = Kpu*f_up/BP; 
                % where pKa is the acidic and pKb is the basic dissocitaion constant. 
           end 
    case 'Arundel'
end
end