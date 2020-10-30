%%
% infusion for IVGTT
if and(strfind(pwd,'IVGTT'),tout(i)<15)
    infRate = 203.37899/15*5;%mmol per 5mn
    harvey = addExchangeRxn(harvey,{'glc_D[bc]'},-1e+6,0);
    harvey = changeRxnBounds(harvey,'EX_glc_D[bc]',-infRate,'b');
end
    %%
    %Lung
    %look for lung reactions
    %findHarveyReactions('Lung',harvey); 
    %%
    organ='Lung';
    % Concentrations
    % Lung_glc_D(c)0
    dglcPulc = ((0-P_10850) ...
            +((P_6650*P_10540*P_6035)-(P_9189*P_10539*P_8543)))/1000;%y(804)
    % Lung_glc_D(e)
    dglcPule = ((((1-P_4463)*P_4464*P_12149*P_6043*P_8554)-((1-P_4460)*(IIf(P_2735,0,P_4465))*P_6043*P_6035)) ...
            +(0-(P_1821*P_6035)) ...
            +(P_2754*P_9170*P_4837*(P_6018-(P_6035/P_4833))) ...
            +(P_2754*(((P_1819+(P_1818*P_1821))*(1-P_9172)*P_6018)+((P_9169*P_4837)*(P_6018-(P_6035/P_4833))*(P_10532/((exp(P_10532))+(-1))))+((0+(((1-P_1818)*P_1821)-P_1819))*(1-P_9173)*P_6018)+((P_9171*P_4837)*(P_6018-(P_6035/P_4833))*(P_10531/((exp(P_10531))+(-1)))))) ...
            +(0-((P_6650*P_10540*P_6035)-(P_9189*P_10539*P_8543))))/1000;%y(800)
    % Lung_g6p
    dg6pPulc = P_10850/1000;%y(806)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcPulc*5;% glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcPule*5;% glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets)))   = dg6pPulc*5;% g6p_D(c)
        
    %Transport    
    met='glc_D[c]';
    oldNumVec = [24871 26532 67180];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_ + P_
    derivCons = ((P_6650*P_10540*P_6035)-(P_9189*P_10539*P_8543))/1000;%y
    %subject constraint
    harvey.b(end) = derivCons*5;

    %Hexokinase  
    met='g6p[c]';
    oldNumVec = [25404 25405 26691];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10850/1000*5;
    
    %print constraints
    dglcPulc = dglcPulc*5
    dglcPule = dglcPule*5
    dg6pPulc = dg6pPulc*5
    transPul = derivCons*5
    hexPul   = P_10850/1000*5
    %%
    %Brain
    %look for brain reactions
    %findHarveyReactions('Brain',harvey); 
%%
    organ='Brain';
    % Concentrations
    % brain_glc_D(c)
    dglcBrainc = ((0-P_10738) ...
            +((P_6250*P_10407*P_5413)-(P_8855*P_10408*P_7981)) ...
            +(P_5396*P_157*P_179*P_5116*P_3841*(((P_5395-P_7981)/P_2771)/(P_181*((1+(P_5395/P_2771))+(1+(P_7981/P_2771))+(-1))))))/1000;%y(75)
    % brain_glc_D(e)
    dglcBraine = ((P_2754*P_8834*P_4684*(P_5395-(P_5413/P_4680))) ...
            +(P_2754*(((P_149+(P_150*P_154))*(1-P_8836)*P_5395)+((P_8833*P_4684)*(P_5395-(P_5413/P_4680))*(P_10399/((exp(P_10399))+(-1))))+((0+(((1-P_150)*P_154)-P_149))*(1-P_8837)*P_5395)+((P_8835*P_4684)*(P_5395-(P_5413/P_4680))*(P_10400/((exp(P_10400))+(-1)))))) ...
            +(0-(P_154*P_5413)) ...
            +(0-((P_6250*P_10407*P_5413)-(P_8855*P_10408*P_7981))) ...
            +(((1-P_3837)*P_3835*P_12113*P_5421*P_7992)-((1-P_3833)*(IIf(P_2735,0,P_3836))*P_5421*P_5413)))/1000;%y(71)
    % brain_g6p
    dg6pBrainc = P_10738/1000;%y(77)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcBrainc*5;% brain_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcBraine*5;% brain_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets)))= dg6pBrainc*5;% brain_g6p

    %Transport
    met='glc_D[c]';
    oldNumVec = [21714 21715 23212 23635]; 
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_ + P_
    derivCons = ((P_6250*P_10407*P_5413)-(P_8855*P_10408*P_7981))/1000;
    harvey.b(end) = derivCons*5;

    %Hexokinase
    met='g6p[c]';
    oldNumVec = [22269 22270 23350];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10738/1000*5;

    %Glut 3 (in Harvey that was SGLT)
    derivCons3 = (P_5396*P_157*P_179*P_5116*P_3841*(((P_5395-P_7981)/P_2771)/(P_181*((1+(P_5395/P_2771))+(1+(P_7981/P_2771))+(-1)))))/1000;
    harvey = changeRxnBounds(harvey,harvey.rxns(23115),derivCons3*5,'b');%GLUT3
    
    %print constraints
    dglcBrainc = dglcBrainc*5
    dglcBraine = dglcBraine*5
    dg6pBrainc = dg6pBrainc*5
    transBrain = derivCons*5
    hexBrain   = P_10738/1000*5
    glut3Brain = derivCons3*5
    %%
    %Heart 
    %look for heart reactions
    %findHarveyReactions('Heart',harvey);
    %%
    organ='Heart';
    %Concentrations
    % Heart_glc_D(c)
    dglcHeartc = ((0-P_10743) ...
            +((P_6298*P_10455*P_5499)-(P_8949*P_10456*P_8095)))/1000;%y(195)
    % heart_glc_D(e)
    dglcHearte = ((P_2754*P_8931*P_4705*(P_5482-(P_5499/P_4701))) ...
            +(P_2754*(((P_416+(P_417*P_418))*(1-P_8929)*P_5482)+((P_8930*P_4705)*(P_5482-(P_5499/P_4701))*(P_10447/((exp(P_10447))+(-1))))+((0+(((1-P_417)*P_418)-P_416))*(1-P_8933)*P_5482)+((P_8932*P_4705)*(P_5482-(P_5499/P_4701))*(P_10448/((exp(P_10448))+(-1)))))) ...
            +(0-((P_6298*P_10455*P_5499)-(P_8949*P_10456*P_8095))) ...
            +(((1-P_3914)*P_3912*P_12125*P_5507*P_8106)-((1-P_3911)*(IIf(P_2735,0,P_3913))*P_5507*P_5499)) ...
            +(0-(P_418*P_5499)))/1000;%y(191)
    % heart_g6p
    dg6pHeartc = P_10743/1000;%y(197)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcHeartc*5;% brain_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcHearte*5;% brain_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets)))   = dg6pHeartc*5;% brain_g6p

    %Transport
    met='glc_D[c]';
    oldNumVec = [35125 35126 36615 36705 37127];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_ + P_
    derivCons = ((P_6298*P_10455*P_5499)-(P_8949*P_10456*P_8095))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;

    %Hexokinase
    met='g6p[c]';
    oldNumVec = [35634 35635 36951];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10743/1000*5;
    
    %print constraints
    dglcHeartc = dglcHeartc*5
    dglcHearte = dglcHearte*5
    dg6pHeartc = dg6pHeartc*5
    transHeart = derivCons*5
    hexHeart   = P_10743/1000*5
    %%
    %Stomach
    %look for Stomach reactions
    %findHarveyReactions('Stomach',harvey);
    %%
    organ='Stomach';
    %Constraints
    %Concentrations
    % stomach_glc_D(c)
    dglcStomachc = ((0-P_10760) ...
            +((P_6450*P_10492*P_5964)-(P_9019*P_10491*P_8351)))/1000;%y(380)
    % stomach_glc_D(e)
    dglcStomache = ((((1-P_4160)*P_4158*P_12133*P_5972*P_8362)-((1-P_4157)*(IIf(P_2735,0,P_4159))*P_5972*P_5964)) ...
            +(0-((P_6450*P_10492*P_5964)-(P_9019*P_10491*P_8351))) ...
            +(0-(P_825*P_5964)) ...
            +(P_2754*P_9032*P_4767*(P_5947-(P_5964/P_4763))) ...
            +(P_2754*(((P_822+(P_820*P_825))*(1-P_9029)*P_5947)+((P_9031*P_4767)*(P_5947-(P_5964/P_4763))*(P_10499/((exp(P_10499))+(-1))))+((0+(((1-P_820)*P_825)-P_822))*(1-P_9030)*P_5947)+((P_9033*P_4767)*(P_5947-(P_5964/P_4763))*(P_10500/((exp(P_10500))+(-1)))))))/1000;%y(376)
    % stomach_g6p
    dg6pStomachc = P_10760/1000;%y(382)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcStomachc*5;% Stomach_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcStomache*5;% Stomach_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets))) = dg6pStomachc*5;% Stomach_g6p

    %Transport
    met='glc_D[c]';
    oldNumVec = [45448 45449 46177];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_ + P_
    derivCons = ((P_6450*P_10492*P_5964)-(P_9019*P_10491*P_8351))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;

    %Hexokinase
    met       ='g6p[c]';
    oldNumVec = [45883 45884 46272];
    harvey    = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10760/1000*5;
    
    %print constraints
    dglcStomachc = dglcStomachc*5
    dglcStomache = dglcStomache*5
    dg6pStomachc = dg6pStomachc*5
    transStomach = derivCons*5
    hexStomach   = P_10760/1000*5
    %%
    %Skin
    %look for Skin reactions
    %findHarveyReactions('Skin',harvey);
    %%
    organ='Skin';
    %Constraints
    %Concentrations
    % skin_glc_D(c)
    dglcSkinc = ((0-P_10855) ...
            +((P_6754*P_10583*P_6145)-(P_9277*P_10584*P_8667)))/1000;%y(939)
    % skin_glc_D(e)
    dglcSkine = ((0-(P_2118*P_6145)) ...
            +(0-((P_6754*P_10583*P_6145)-(P_9277*P_10584*P_8667))) ...
            +(((1-P_4549)*P_4551*P_12161*P_6153*P_8678)-((1-P_4547)*(IIf(P_2735,0,P_4552))*P_6153*P_6145)) ...
            +(P_2754*P_9289*P_4884*(P_6128-(P_6145/P_4880))) ...
            +(P_2754*(((P_2120+(P_2119*P_2118))*(1-P_9291)*P_6128)+((P_9293*P_4884)*(P_6128-(P_6145/P_4880))*(P_10591/((exp(P_10591))+(-1))))+((0+(((1-P_2119)*P_2118)-P_2120))*(1-P_9292)*P_6128)+((P_9290*P_4884)*(P_6128-(P_6145/P_4880))*(P_10592/((exp(P_10592))+(-1)))))))/1000;%y(935)
    % skin_g6p
    dg6pSkinc = P_10855/1000;%y(941)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcSkinc*5;% skin_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcSkine*5;% skin_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets))) = dg6pSkinc*5;% skin_g6p_D(c)
    
    %Transport
    met='glc_D[c]';
    oldNumVec = [30179 30180 31195];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_7392 + P_7393
    derivCons = ((P_6754*P_10583*P_6145)-(P_9277*P_10584*P_8667))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;
    
    %Hexokinase
    met='g6p[c]';
    oldNumVec = [30599 30600 31306];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10855/1000*5;
    
    %print constraints
    dglcSkinc = dglcSkinc*5
    dglcSkine = dglcSkine*5
    dg6pSkinc = dg6pSkinc*5
    skinTrans = derivCons*5
    hexSkin   = P_10855/1000*5
    %%
    %Gonads (Testis)
    %look for Testis reactions
    %findHarveyReactions('Testis',harvey);
    %%
    organ='Testis';
    %Constraints
    %Concentrations
    % Testis_glc_D(c)
    dglcTestisc = ((0-P_10742) ...
            +((P_6274*P_10431*P_5471)-(P_8899*P_10432*P_8061)))/1000;%y(159)
    % Testis_glc_D(e)
    dglcTestise = ((((1-P_3887)*P_3889*P_12121*P_5479*P_8072)-((1-P_3884)*(IIf(P_2735,0,P_3888))*P_5479*P_5471)) ...
            +(0-((P_6274*P_10431*P_5471)-(P_8899*P_10432*P_8061))) ...
            +(0-(P_345*P_5471)) ...
            +(P_2754*P_8912*P_4699*(P_5454-(P_5471/P_4695))) ...
            +(P_2754*(((P_344+(P_346*P_345))*(1-P_8909)*P_5454)+((P_8911*P_4699)*(P_5454-(P_5471/P_4695))*(P_10440/((exp(P_10440))+(-1))))+((0+(((1-P_346)*P_345)-P_344))*(1-P_8910)*P_5454)+((P_8913*P_4699)*(P_5454-(P_5471/P_4695))*(P_10439/((exp(P_10439))+(-1)))))))/1000;%y(155)
    % Testis_g6p
    dg6pTestisc = P_10742/1000;%y(161)
    
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcTestisc*5;% 
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcTestise*5;% 
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets))) = dg6pTestisc*5;% 
    
    %Transport
    met='glc_D[c]';
    oldNumVec = [47141 47142 48207 48310 48676];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_7278 + P_7279
    derivCons = ((P_6274*P_10431*P_5471)-(P_8899*P_10432*P_8061))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;
    
    %Hexokinase
    met='g6p[c]';
    oldNumVec = [47605 47606 48436];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10742/1000*5;
    
    %print constraints
     dglcTestisc = dglcTestisc*5
     dglcTestise = dglcTestise*5
     dg6pTestisc = dg6pTestisc*5
     testisTrans = derivCons*5
     hexTestis   = P_10742/1000*5
    %%
    %Spleen
    %look for Spleen reactions
    %findHarveyReactions('Spleen',harvey);
    %%
    organ='Spleen';
    %Constraints
    %Concentrations
    % spleen_glc_D(c)
    dglcSpleenc = ((0-P_10856) ...
            +((P_6770*P_10599*P_6173)-(P_9311*P_10600*P_8701)))/1000;%y(975)
    % spleen_glc_D(e)
    dglcSpleene =             ((((1-P_4578)*P_4579*P_12165*P_6181*P_8712)-((1-P_4574)*(IIf(P_2735,0,P_4577))*P_6181*P_6173)) ...
            +(0-((P_6770*P_10599*P_6173)-(P_9311*P_10600*P_8701))) ...
            +(0-(P_2189*P_6173)) ...
            +(P_2754*P_9323*P_4891*(P_6156-(P_6173/P_4887))) ...
            +(P_2754*(((P_2191+(P_2192*P_2189))*(1-P_9324)*P_6156)+((P_9322*P_4891)*(P_6156-(P_6173/P_4887))*(P_10607/((exp(P_10607))+(-1))))+((0+(((1-P_2192)*P_2189)-P_2191))*(1-P_9325)*P_6156)+((P_9321*P_4891)*(P_6156-(P_6173/P_4887))*(P_10608/((exp(P_10608))+(-1)))))))/1000;%y(971)
    % spleen_g6p 
    dg6pSpleenc = P_10856/1000;%y(977)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcSpleenc*5;% spleen_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcSpleene*5;% spleen_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets))) = dg6pSpleenc*5;% spleen_g6p

    %Transport
    met='glc_D[c]';
    oldNumVec = [20117];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_7404 + P_7405
    derivCons = ((P_6770*P_10599*P_6173)-(P_9311*P_10600*P_8701))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;

    %Hexokinase
    met='g6p[c]';
    oldNumVec = [20473 20474 20830];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10856/1000*5;
    
    %print constraints
    dglcSpleenc = dglcSpleenc*5
    dglcSpleene = dglcSpleene*5
    dg6pSpleenc = dg6pSpleenc*5
    transHex    = derivCons*5
    hexSpleen   = P_10856/1000*5
    %%
%     %Kidney
%     look for Kidney reactions
%     findHarveyReactions('Kidney',harvey);
%     check if glucose in urine exists
%     newMet  = 'glc_D[u]';
%     [o1,r1] = ismember(newMet,harvey.mets);
%     c = find(harvey.S(r1,:));
%     for i = 1:length(c)
%         c(i)
%         printRxnFormula(harvey,harvey.rxns([c(i)]));
%     end
%%
    %set glucose secretion in urine for diabetic
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'glc_D[u]')),harvey.mets))) = (P_2754*P_5509*P_2762*P_519*P_5511)/1000*5;%y(253)
    glcUrine = (P_2754*P_5509*P_2762*P_519*P_5511)/1000*5
    %%
    organ='Kidney';
    %Constraints
    %Concentrations
    % Kidney_glc_D(c) 
    dglcKidneyc = ((0-P_10744) ...
            +((P_6314*P_10475*P_5529)-(P_8983*P_10476*P_8129)))/1000;%y(231)
    % Kidney_glc_D(e)
    dglcKidneye = ((0-(P_488*P_5529)) ...
            +(P_2754*P_8963*P_4712*(P_5511-(P_5529/P_4708))) ...
            +(P_2754*(((P_491+(P_490*P_488))*(1-P_8965)*P_5511)+((P_8962*P_4712)*(P_5511-(P_5529/P_4708))*(P_10467/((exp(P_10467))+(-1))))+((0+(((1-P_490)*P_488)-P_491))*(1-P_8961)*P_5511)+((P_8964*P_4712)*(P_5511-(P_5529/P_4708))*(P_10468/((exp(P_10468))+(-1)))))) ...
            +(((1-P_3939)*P_3942*P_12129*P_5537*P_8140)-((1-P_3938)*(IIf(P_2735,0,P_3943))*P_5537*P_5529)) ...
            +(0-((P_6314*P_10475*P_5529)-(P_8983*P_10476*P_8129))))/1000;%y(227)
    % Kidney_g6p
    dg6pKidneyc = P_10744/1000;%y(233)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcKidneyc*5;% Kidney_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcKidneye*5;% Kidney_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets))) = dg6pKidneyc*5;% Kidney_g6p

    %Transport
    met='glc_D[c]';
    oldNumVec = [15607 15608 17354 17491 18108];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_7296 + P_7297
    derivCons = ((P_6314*P_10475*P_5529)-(P_8983*P_10476*P_8129))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;

    %Hexokinase
    met='g6p[c]';
    oldNumVec = [16264 16263 17673];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10744/1000*5;
    
    %print constraints
    dglcKidneyc = dglcKidneyc*5
    dglcKidneye = dglcKidneye*5
    dg6pKidneyc = dg6pKidneyc*5
    transKidney = derivCons*5
    KidneyHex   = P_10744/1000*5
    %%
    %Gallbladder (to be added at further stages)
    %look for Gallbladder reactions
    %findHarveyReactions('Gall',harvey);
    %%
    %Adipocytes
    %look for Adipocyte reactions
    %findHarveyReactions('Adipocytes',harvey);
    %%
    organ='Adipocytes';
    %Constraints
    %Concentrations
    % Adipocyte_glc_D(c)
    dglcAdipocytec = ((0-P_10739) ...
            +(P_2775*(y(108)^(2*P_2294))*P_5117*P_3874*(((P_6258*0)+P_5441)/(P_276*(P_2774+(P_6258*0)+P_5441)))) ...
            +((P_6258*P_10415*P_5441)-(P_8866*P_10416*P_8015)))/1000;%y(112)
    % Adipocyte_glc_D(e)
    dglcAdipocytee = ((0-(P_230*P_5441)) ...
            +(0-(P_2775*(y(108)^(2*P_2294))*P_5117*P_3874*(((P_6258*0)+P_5441)/(P_276*(P_2774+(P_6258*0)+P_5441))))) ...
            +(0-((P_6258*P_10415*P_5441)-(P_8866*P_10416*P_8015))) ...
            +(((1-P_3861)*P_3863*P_12117*P_5451*P_8038)-((1-P_3859)*(IIf(P_2735,0,P_3860))*P_5451*P_5441)) ...
            +(P_2754*P_8880*P_4692*(P_5424-(P_5441/P_4688))) ...
            +(P_2754*(((P_234+(P_231*P_230))*(1-P_8877)*P_5424)+((P_8879*P_4692)*(P_5424-(P_5441/P_4688))*(P_10423/((exp(P_10423))+(-1))))+((0+(((1-P_231)*P_230)-P_234))*(1-P_8878)*P_5424)+((P_8881*P_4692)*(P_5424-(P_5441/P_4688))*(P_10424/((exp(P_10424))+(-1)))))))/1000;%y(107)
    % Adipocyte_g6p
    dg6pAdipocytec = P_10739/1000;%y(117)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcAdipocytec*5;% Adipocyte_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcAdipocytee*5;% Adipocyte_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets))) = dg6pAdipocytec*5;% Adipocyte_g6p

    %Transport
    met='glc_D[c]';
    oldNumVec = [32110];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_7260 + P_7261
    derivCons = ((P_6258*P_10415*P_5441)-(P_8866*P_10416*P_8015))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;

    %Glut 4 (not exactly glu 4 in harvey)
    derivCons3 = (P_2775*(y(108)^(2*P_2294))*P_5117*P_3874*(((P_6258*0)+P_5441)/(P_276*(P_2774+(P_6258*0)+P_5441))))/1000;
    harvey = changeRxnBounds(harvey,harvey.rxns(32697),derivCons3*5,'b');%GLUT4

    %Hexokinase
    met='g6p[c]';
    oldNumVec = [32822 32442 32441];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10739/1000*5;
    
    %print constraints
    dglcAdipocytec = dglcAdipocytec*5
    dglcAdipocytee = dglcAdipocytee*5
    dg6pAdipocytec = dg6pAdipocytec*5
    transAdipocyte = derivCons*5
    glutAdipocyte  = derivCons3*5
    hexAdipocyte   = P_10739/1000*5
    %%
    %Muscle
    %look for Muscle reactions
    %findHarveyReactions('Muscle',harvey);
    %%
    organ='Muscle';
    %Constraints
    %Concentrations
    % Muscle_glc_D(c)
    dglcMusclec = ((0-P_10851) ...
            +(P_2775*(y(837)^(2*P_2294))*P_5117*P_4501*(((P_6658*0)+P_6063)/(P_1935*(P_2774+(P_6658*0)+P_6063)))) ...
            +((P_6658*P_10547*P_6063)-(P_9203*P_10548*P_8577)))/1000;%y(841)
    % Muscle_glc_D(e)
    dglcMusclee = ((((1-P_4489)*P_4487*P_12153*P_6073*P_8600)-((1-P_4486)*(IIf(P_2735,0,P_4488))*P_6073*P_6063)) ...
            +(0-(P_2775*(y(837)^(2*P_2294))*P_5117*P_4501*(((P_6658*0)+P_6063)/(P_1935*(P_2774+(P_6658*0)+P_6063))))) ...
            +(0-((P_6658*P_10547*P_6063)-(P_9203*P_10548*P_8577))) ...
            +(0-(P_1895*P_6063)) ...
            +(P_2754*P_9215*P_4845*(P_6046-(P_6063/P_4841))) ...
            +(P_2754*(((P_1892+(P_1891*P_1895))*(1-P_9213)*P_6046)+((P_9214*P_4845)*(P_6046-(P_6063/P_4841))*(P_10555/((exp(P_10555))+(-1))))+((0+(((1-P_1891)*P_1895)-P_1892))*(1-P_9217)*P_6046)+((P_9216*P_4845)*(P_6046-(P_6063/P_4841))*(P_10556/((exp(P_10556))+(-1)))))))/1000;%y(836)
    % Muscle_g6p
    dg6pMusclec = P_10851/1000;%y(846)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcMusclec*5;% Muscle_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcMusclee*5;% Muscle_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets))) = dg6pMusclec*5;% Muscle_g6p

    %Transport
    met='glc_D[c]';
    oldNumVec = [27814 27815];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_7365 + P_7366
    derivCons = ((P_6658*P_10547*P_6063)-(P_9203*P_10548*P_8577))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;

    %Glut 4
    derivCons3 = (P_2775*(y(837)^(2*P_2294))*P_5117*P_4501*(((P_6658*0)+P_6063)/(P_1935*(P_2774+(P_6658*0)+P_6063))))/1000;
    harvey = changeRxnBounds(harvey,harvey.rxns(29084),derivCons3*5,'b');%GLUT4

    %Hexokinase
    met='g6p[c]';
    oldNumVec = [29215 28298 28299];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10851/1000*5;
    
    %print constraints
    dglcMusclec = dglcMusclec*5
    dglcMusclee = dglcMusclee*5
    dg6pMusclec = dg6pMusclec*5
    transMuscle = derivCons*5
    glut4Muscle = derivCons3*5
    hexMuscle   = P_10851/1000*5
    %%
    %Liver
    %look for Liver reactions
    %findHarveyReactions('Liver',harvey);
    %%
    organ='Liver';
    %Constraints
    %Concentrations
    % Liver_glc_D(c) 
    dglcLiverc = ((0-(P_10845*y(764)*P_10257*P_8491)) ...
            +(P_12108*y(762)*P_8489*P_10256) ...
            +(0-(P_1773*P_5115*P_6013*(((P_8493/P_2767)-(P_6003/P_2768))/(P_1764*((1+(P_8493/P_2767))+(1+(P_6003/P_2768))+(-1)))))) ...
            +((P_6634*P_10523*P_6003)-(P_9159*P_10524*P_8493)))/1000;%y(754)
    % Liver_glc_D(e)
    dglcLivere = ((0-(P_1673*P_6003)) ...
            +(P_2754*P_9138*P_4830*(P_5986-(P_6003/P_4826))) ...
            +(P_2754*(((P_1670+(P_1669*P_1673))*(1-P_9140)*P_5986)+((P_9137*P_4830)*(P_5986-(P_6003/P_4826))*(P_10515/((exp(P_10515))+(-1))))+((0+(((1-P_1669)*P_1673)-P_1670))*(1-P_9141)*P_5986)+((P_9139*P_4830)*(P_5986-(P_6003/P_4826))*(P_10516/((exp(P_10516))+(-1)))))) ...
            +(((1-P_4437)*P_4435*P_12145*P_6014*P_8520)-((1-P_4434)*(IIf(P_2735,0,P_4436))*P_6014*P_6003)) ...
            +(P_1773*P_5115*P_6013*(((P_8493/P_2767)-(P_6003/P_2768))/(P_1764*((1+(P_8493/P_2767))+(1+(P_6003/P_2768))+(-1))))) ...
            +(0-((P_6634*P_10523*P_6003)-(P_9159*P_10524*P_8493))))/1000;%y(748)
    % Liver_glycogen
    dglyLiverc = ((P_10845*y(764)*P_10257*P_8491) ...
            +(0-(P_12108*y(762)*P_8489*P_10256)))/1000;%y(761)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcLiverc*5;% Liver_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcLivere*5;% Liver_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glygn2[c]'])),harvey.mets))) = dglyLiverc*5;% Liver_ggn (Glycogen)

    %Glut 2
    derivCons3 = (P_1773*P_5115*P_6013*(((P_8493/P_2767)-(P_6003/P_2768))/(P_1764*((1+(P_8493/P_2767))+(1+(P_6003/P_2768))+(-1)))))/1000;
    harvey = changeRxnBounds(harvey,harvey.rxns(8740),derivCons3*5,'b');%GLUT2
    
    %Transport
    met='glc_D[c]';
    oldNumVec = [6445 6446];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_7342 + P_7341
    derivCons = ((P_6634*P_10523*P_6003)-(P_9159*P_10524*P_8493))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;
    
    %Glycogen phosphorylase
    derivCons4 = (P_12108*y(762)*P_8489*P_10256)/1000;
    harvey = changeRxnBounds(harvey,harvey.rxns(6459),derivCons4*5,'b');
    
    %Hexokinase // Glycogen synthase
    met='g6p[c]';
    oldNumVec = [7169 7170 9072];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = (P_10845*y(764)*P_10257*P_8491)/1000*5;
    
    %print constraints
    dglcLiverc = dglcLiverc*5
    dglcLivere = dglcLivere*5
    dglyLiverc = dglyLiverc*5
    glut2Liver = derivCons3*5
    transLiver = derivCons*5
    glycoLiver = derivCons4*5
    hexLiver = (P_10845*y(764)*P_10257*P_8491)/1000*5
    %%
    %Pancreas
    %look for Pancreas reactions
    %findHarveyReactions('Pancreas',harvey);
    %%
    organ='Pancreas';
    %Constraints
    %Concentrations
    % Pancreas_glc_D(c)
    dglcPancreasc = ((0-P_10854) ...
            +((P_6682*P_10571*P_6093)-(P_9253*P_10572*P_8629)))/1000;%y(893)
    % Pancreas_glc_D(e)
    dglcPancrease = ((P_2754*P_9234*P_4850*(P_6076-(P_6093/P_4846))) ...
            +(P_2754*(((P_2003+(P_2002*P_2001))*(1-P_9237)*P_6076)+((P_9233*P_4850)*(P_6076-(P_6093/P_4846))*(P_10563/((exp(P_10563))+(-1))))+((0+(((1-P_2002)*P_2001)-P_2003))*(1-P_9236)*P_6076)+((P_9235*P_4850)*(P_6076-(P_6093/P_4846))*(P_10564/((exp(P_10564))+(-1)))))) ...
            +(0-(P_2001*P_6093)) ...
            +(0-((P_6682*P_10571*P_6093)-(P_9253*P_10572*P_8629))) ...
            +(((1-P_4514)*P_4512*P_12157*P_6108*P_8640)-((1-P_4510)*(IIf(P_2735,0,P_4513))*P_6108*P_6093)))/1000;%y(884)
    % Pancreas_g6p 
    dg6pPancreasc = P_10854/1000;%y(895)
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[c]'])),harvey.mets))) = dglcPancreasc*5;% Pancreas_glc_D(c)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_glc_D[e]'])),harvey.mets))) = dglcPancrease*5;% Pancreas_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,[organ '_g6p[c]'])),harvey.mets))) = dg6pPancreasc*5;% Pancreas_g6p
    
    %Transport
    met='glc_D[c]';
    oldNumVec = [12975 66010];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    %compute cosntraint = P_7383 + P_7384
    derivCons = ((P_6682*P_10571*P_6093)-(P_9253*P_10572*P_8629))/1000;
    %subject constraint
    harvey.b(end) = derivCons*5;

    %Hexokinase
    met='g6p[c]';
    oldNumVec = [13439 13440 13844];
    harvey = addMetabolite2Harvey(harvey,oldNumVec,organ,met);
    harvey.b(end) = P_10854/1000*5;
    
    %print constraints
    dglcPancreasc = dglcPancreasc*5
    dglcPancrease = dglcPancrease*5
    dg6pPancreasc = dg6pPancreasc*5
    transPancreas = derivCons*5
    hexPancreas   = P_10854/1000*5
    %%
%     Portal Vein
%     look for Portal vein reactions
%     [o2,r2] = ismember('glc_D[bp]',harvey.mets);
%     sprintf('glc_D[bp]')
%     a = find(harvey.S(r2,:));
%     for i = 1:length(a)
%         a(i)
%         affich = printRxnFormula(harvey,harvey.rxns([a(i)]));
%     end
    %%
    %add reactions (spleen, stomach, and pancreas are connected to the portal vein in GIM but not in Harvey)
    %Spleen
    harvey = addReaction(harvey,'Spleen_BL_secretion1',{'glc_D[bp]','Spleen_glc_D[bpS]'},[1,-1],false,[],[],[],[],'',[],[],[]);
    harvey = addReaction(harvey,'Spleen_BL_secretion2',{'Spleen_glc_D[c]','Spleen_glc_D[bpS]'},[-1,1],true,[],[],[],[],'',[],[],[]);
    %Stomach
    harvey = addReaction(harvey,'Stomach_BL_secretion1',{'glc_D[bp]','Stomach_glc_D[bpSt]'},[1,-1],false,[],[],[],[],'',[],[],[]);
    harvey = addReaction(harvey,'Stomach_BL_secretion2',{'Stomach_glc_D[c]','Stomach_glc_D[bpSt]'},[-1,1],true,[],[],[],[],'',[],[],[]);
    %Pancreas
    harvey = addReaction(harvey,'Pancreas_BL_secretion1',{'glc_D[bp]','Pancreas_glc_D[bpP]'},[1,-1],false,[],[],[],[],'',[],[],[]);
    harvey = addReaction(harvey,'Pancreas_BL_secretion2',{'Pancreas_glc_D[c]','Pancreas_glc_D[bpP]'},[-1,1],true,[],[],[],[],'',[],[],[]);
    %%
    %Constraints
    %Concentrations
    % PortalVein_glc_D(e)
    dglcPortalVeine = (P_11719 ...
            +P_11735 ...
            +(((P_4164*(1-P_40))-P_825)*P_5947) ...
            +P_11751 ...
            +((((1-P_893)*P_4194*(1-P_40))-P_11064)*P_11552) ...
            +P_11779 ...
            +((((1-P_1287)*P_4310*(1-P_40))-P_11232)*P_11629) ...
            +(((P_4517*(1-P_40))-P_2001)*P_6076) ...
            +P_11815 ...
            +(0-(P_4877*P_2754*((P_9835*P_6111)-(P_9833*(P_6118/P_7396))))) ...
            +(0-(((P_6109*(1-P_40))-P_11708)*P_6111)) ...
            +(((P_4581*(1-P_40))-P_2189)*P_6156) ...
            +P_11824 ...
            +P_11844 ...
            +P_11876 ...
            +P_11888 ...
            +P_11892 ...
            +P_11904)/1000;%y(915)
    % Spleen_portalblood_glc_D(bpS)
    dglcPortalbloodbpS = ((0-(P_4886*P_2754*((P_9835*P_6156)-(P_9833*(P_6164/P_7396))))) ...
            +(0-((P_4574*(IIf(P_2735,0,P_4577))*P_6181*P_6156)-(P_4578*P_4579*P_12349*P_6181*P_8712))) ...
            +(0-(P_2754*P_9323*P_4891*(P_6156-(P_6173/P_4887)))) ...
            +(0-(P_2754*(((P_2191+(P_2192*P_2189))*(1-P_9324)*P_6156)+((P_9322*P_4891)*(P_6156-(P_6173/P_4887))*(P_10607/((exp(P_10607))+(-1))))+((0+(((1-P_2192)*P_2189)-P_2191))*(1-P_9325)*P_6156)+((P_9321*P_4891)*(P_6156-(P_6173/P_4887))*(P_10608/((exp(P_10608))+(-1))))))) ...
            +(0-(((P_4581*(1-P_40))-P_2189)*P_6156)) ...
            +(P_4581*(1-P_40)*P_5350))/1000;%y(961)
    % SI_portalblood_glc_D(bpI)  
    DuodenumAmount1     = yout(i,438)*P_10769;
    UpperJejunumAmount1 = yout(i,462)*P_10776;
    LowerJejunumAmount1 = yout(i,486)*P_10783;
    UpperIleumAmount1   = yout(i,510)*P_10790;
    LowerIleumAmount1   = yout(i,534)*P_10797;
    SerosaAmount1       = yout(i,402)*P_11066;
    dglcPortalbloodbpI1 = (DuodenumAmount1 + UpperJejunumAmount1 +...
       LowerJejunumAmount1 + UpperIleumAmount1 + LowerIleumAmount1 +...
       SerosaAmount1) / 1000;
    DuodenumAmount2     = yout(i+1,438)*P_10769;
    UpperJejunumAmount2 = yout(i+1,462)*P_10776;
    LowerJejunumAmount2 = yout(i+1,486)*P_10783;
    UpperIleumAmount2   = yout(i+1,510)*P_10790;
    LowerIleumAmount2   = yout(i+1,534)*P_10797;
    SerosaAmount2       = yout(i+1,402)*P_11066;
    dglcPortalbloodbpI2 = (DuodenumAmount2 + UpperJejunumAmount2 +...
       LowerJejunumAmount2 + UpperIleumAmount2 + LowerIleumAmount2 +...
       SerosaAmount2) / 1000;
    dglcPortalbloodbpI = dglcPortalbloodbpI2 - dglcPortalbloodbpI1;
    % Stomach_portalblood_glc_D(bpSt)
    dglcPortalbloodbpSt = ((0-(((P_4164*(1-P_40))-P_825)*P_5947)) ...
            +(0-((P_4157*(IIf(P_2735,0,P_4159))*P_5972*P_5947)-(P_4160*P_4158*P_12285*P_5972*P_8362))) ...
            +(0-(P_2754*P_9032*P_4767*(P_5947-(P_5964/P_4763)))) ...
            +(0-(P_2754*(((P_822+(P_820*P_825))*(1-P_9029)*P_5947)+((P_9031*P_4767)*(P_5947-(P_5964/P_4763))*(P_10499/((exp(P_10499))+(-1))))+((0+(((1-P_820)*P_825)-P_822))*(1-P_9030)*P_5947)+((P_9033*P_4767)*(P_5947-(P_5964/P_4763))*(P_10500/((exp(P_10500))+(-1))))))) ...
            +(0-(P_4768*P_2754*((P_9835*P_5947)-(P_9833*(P_5955/P_7396))))) ...
            +(P_4164*(1-P_40)*P_5350))/1000;%y(366)
    % Liver_portalblood_glc_D(bpL)
    dglcPortalbloodbpL = ((0-(P_4825*P_2754*((P_9835*P_5986)-(P_9833*(P_5994/P_7396))))) ...
            +(0-((P_4434*(IIf(P_2735,0,P_4436))*P_6014*P_5986)-(P_4437*P_4435*P_12305*P_6014*P_8520))) ...
            +(0-(P_2754*P_9138*P_4830*(P_5986-(P_6003/P_4826)))) ...
            +(0-(P_2754*(((P_1670+(P_1669*P_1673))*(1-P_9140)*P_5986)+((P_9137*P_4830)*(P_5986-(P_6003/P_4826))*(P_10515/((exp(P_10515))+(-1))))+((0+(((1-P_1669)*P_1673)-P_1670))*(1-P_9141)*P_5986)+((P_9139*P_4830)*(P_5986-(P_6003/P_4826))*(P_10516/((exp(P_10516))+(-1))))))) ...
            +(0-((((P_4441+P_6109)*(1-P_40))-(P_1673+P_11708))*P_5986)) ...
            +(((P_6109*(1-P_40))-P_11708)*P_6111) ...
            +(P_4441*(1-P_40)*P_5350))/1000;%y(738)
    % Pancreas_portalblood_glc_D(bpP)
    dglcPortalbloodbpP = ((0-(((P_4517*(1-P_40))-P_2001)*P_6076)) ...
            +(0-(P_2754*P_9234*P_4850*(P_6076-(P_6093/P_4846)))) ...
            +(0-(P_2754*(((P_2003+(P_2002*P_2001))*(1-P_9237)*P_6076)+((P_9233*P_4850)*(P_6076-(P_6093/P_4846))*(P_10563/((exp(P_10563))+(-1))))+((0+(((1-P_2002)*P_2001)-P_2003))*(1-P_9236)*P_6076)+((P_9235*P_4850)*(P_6076-(P_6093/P_4846))*(P_10564/((exp(P_10564))+(-1))))))) ...
            +(0-((P_4510*(IIf(P_2735,0,P_4513))*P_6108*P_6076)-(P_4514*P_4512*P_12329*P_6108*P_8640))) ...
            +(0-(P_4851*P_2754*((P_9835*P_6076)-(P_9833*(P_6084/P_7396))))) ...
            +(P_4517*(1-P_40)*P_5350))/1000;%y(874)
    % Colon_portalblood_glc_D(bpC)
    cecumAmount1  = yout(i,606)*P_11240;
    colonAAmount1 = yout(i,630)*P_11270;
    colonTAmount1 = yout(i,642)*P_11300;
    colonDAmount1 = yout(i,666)*P_11330;
    colonSAmount1 = yout(i,690)*P_11360;
    SerosaAmount1 = yout(i,558)*P_11628;
    dglcPortalbloodbpC1 = (cecumAmount1 + colonAAmount1 + colonTAmount1 +...
        colonDAmount1 + colonSAmount1 + SerosaAmount1) / 1000;
    cecumAmount2  = yout(i+1,606)*P_11240;
    colonAAmount2 = yout(i+1,630)*P_11270;
    colonTAmount2 = yout(i+1,642)*P_11300;
    colonDAmount2 = yout(i+1,666)*P_11330;
    colonSAmount2 = yout(i+1,690)*P_11360;
    SerosaAmount2 = yout(i+1,558)*P_11628;
    dglcPortalbloodbpC2 = (cecumAmount2 + colonAAmount2 + colonTAmount2 +...
       colonDAmount2 + colonSAmount2 + SerosaAmount2) / 1000;
    dglcPortalbloodbpC = dglcPortalbloodbpC2 - dglcPortalbloodbpC1;
    % set b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'Spleen_glc_D[bpS]')),harvey.mets)))   = dglcPortalbloodbpS*5 ;% Spleen_portalblood_glc_D(bpS)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'Stomach_glc_D[bpSt]')),harvey.mets))) = dglcPortalbloodbpSt*5;% Stomach_portalblood_glc_D(bpSt)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'Pancreas_glc_D[bpP]')),harvey.mets))) = dglcPortalbloodbpP*5 ;% Pancreas_portalblood_glc_D(bpP)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'glc_D[bp]')),harvey.mets)))           = dglcPortalVeine*5    ;% PortalVein_glc_D(e)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'sIEC_glc_D[bpI]')),harvey.mets)))     = dglcPortalbloodbpI ;% SI_portalblood_glc_D(bpI)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'Liver_glc_D[bpL]')),harvey.mets)))    = dglcPortalbloodbpL*5 ;% Liver_portalblood_glc_D(bpL)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'Colon_glc_D[bpC]')),harvey.mets)))    = dglcPortalbloodbpC ;% Colon_portalblood_glc_D(bpC)
    
    %print constraints
    dglcPortalbloodbpS  = dglcPortalbloodbpS*5
    dglcPortalbloodbpSt = dglcPortalbloodbpSt*5
    dglcPortalbloodbpP  = dglcPortalbloodbpP*5
    dglcPortalVeine     = dglcPortalVeine*5
    dglcPortalbloodbpI  
    dglcPortalbloodbpL 
    dglcPortalbloodbpL  = dglcPortalbloodbpL*5
    dglcPortalbloodbpC  
    %%
%     look for sIEC reactions
%     organ = 'sIEC';
%     g6p  = strcat(organ,'_g6p[c]');
%     glce = strcat(organ,'_glc_D[bpI]');
%     glcc =strcat(organ,'_glc_D[c]');
%     [o1,r1] = ismember(g6p,harvey.mets);
%     [o2,r2] = ismember(glce,harvey.mets);
%     [o3,r3] = ismember(glcc,harvey.mets);
% 
%     %glc_D(e)
%     sprintf('glc_D(e)')
%     a = find(harvey.S(r2,:));
%     for i = 1:length(a)
%         a(i)
%         affich = printRxnFormula(harvey,harvey.rxns([a(i)]));
%     end
%     glc_D(c)
%     sprintf('glc_D(c)')
%     b = find(harvey.S(r3,:));
%     for i = 1:length(b)
%         b(i)
%         affich = printRxnFormula(harvey,harvey.rxns([b(i)]));
%     end
%     %Brain_g6p
%     sprintf('g6p(c)')
%     c = find(harvey.S(r1,:));
%     for i = 1:length(c)
%         c(i)
%         printRxnFormula(harvey,harvey.rxns([c(i)]));
%     end
    %%
    %Constraints
    %Concentrations
    % sIEC_g6p
    dg6pSiecc = P_11921/1000;%y(418)
    % sIEC_glc(C)
    DuodenumAmount1     = yout(i,450)*P_10773;
    UpperJejunumAmount1 = yout(i,474)*P_10780;
    LowerJejunumAmount1 = yout(i,498)*P_10787;
    UpperIleumAmount1   = yout(i,522)*P_10794;
    LowerIleumAmount1   = yout(i,546)*P_10801;
    SerosaAmount1       = yout(i,416)*P_11070;
    dglcPortalbloodbpI1 = (DuodenumAmount1 + UpperJejunumAmount1 +...
        LowerJejunumAmount1 + UpperIleumAmount1 + LowerIleumAmount1 +...
        SerosaAmount1) / 1000;
    DuodenumAmount2     = yout(i+1,450)*P_10773;
    UpperJejunumAmount2 = yout(i+1,474)*P_10780;
    LowerJejunumAmount2 = yout(i+1,498)*P_10787;
    UpperIleumAmount2   = yout(i+1,522)*P_10794;
    LowerIleumAmount2   = yout(i+1,546)*P_10801;
    SerosaAmount2       = yout(i+1,416)*P_11070;
    dglcPortalbloodbpI2 = (DuodenumAmount2 + UpperJejunumAmount2 +...
        LowerJejunumAmount2 + UpperIleumAmount2 + LowerIleumAmount2 +...
        SerosaAmount2) / 1000;
    dglcSiecc = dglcPortalbloodbpI2 - dglcPortalbloodbpI1;
    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'sIEC_g6p[c]')),harvey.mets)))   = dg6pSiecc*5;% Siec_g6pc
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'sIEC_glc_D[c]')),harvey.mets))) = dglcSiecc;% Siec_glcc
    
    %print constraints
    dg6pSiecc = dg6pSiecc*5
    dglcSiecc
    %%
    %look for Colon reactions
    %%
    %Colon_g6p
    dg6pColonc = P_11940/1000;%y(574)
    %sIEC_glc(C) 
    cecumAmount1  = yout(i,606)*P_10809;
    colonAAmount1 = yout(i,630)*P_10816;
    colonTAmount1 = yout(i,654)*P_10823;
    colonDAmount1 = yout(i,678)*P_10830;
    colonSAmount1 = yout(i,702)*P_10837;
    SerosaAmount1 = yout(i,572)*P_11238;
    dglcColonC1 = (cecumAmount1 + colonAAmount1 + colonTAmount1 +...
        colonDAmount1 + colonSAmount1 + SerosaAmount1) / 1000;
    cecumAmount2  = yout(i+1,606)*P_10809;
    colonAAmount2 = yout(i+1,630)*P_10816;
    colonTAmount2 = yout(i+1,654)*P_10823;
    colonDAmount2 = yout(i+1,678)*P_10830;
    colonSAmount2 = yout(i+1,702)*P_10837;
    SerosaAmount2 = yout(i+1,572)*P_11238;
    dglcColonC2 = (cecumAmount2 + colonAAmount2 + colonTAmount2 +...
        colonDAmount2 + colonSAmount2 + SerosaAmount2) / 1000;
    dglcColonC = dglcColonC2 - dglcColonC1;
    %fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'Colon_g6p[c]')),harvey.mets)))   = dg6pColonc*5;% Colon_g6p
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'Colon_glc_D[c]')),harvey.mets))) = dglcColonC;%finite sum between 2 pts -> in 5mn
    
    %print constraints
    dglcColonC = dglcColonC
    dg6pColonc = dg6pColonc*5
    %%
    %look for Lumen reactions
    %%
    %Constraints
    %Concentrations
    % Lumen_glc_D(c)

    DuodenumAmount1     = yout(i,263)*P_5562;
    UpperJejunumAmount1 = yout(i,272)*P_5597;
    LowerJejunumAmount1 = yout(i,281)*P_5632;
    UpperIleumAmount1   = yout(i,290)*P_5667;
    LowerIleumAmount1   = yout(i,299)*P_5700;
    dglcLumenc1 = (DuodenumAmount1 + UpperJejunumAmount1 +...
    LowerJejunumAmount1 + UpperIleumAmount1 + LowerIleumAmount1 ...
    ) / 1000;
    DuodenumAmount2     = yout(i+1,263)*P_5562;
    UpperJejunumAmount2 = yout(i+1,272)*P_5597;
    LowerJejunumAmount2 = yout(i+1,281)*P_5632;
    UpperIleumAmount2   = yout(i+1,290)*P_5667;
    LowerIleumAmount2   = yout(i+1,299)*P_5700;
    dglcLumenc2 = (DuodenumAmount2 + UpperJejunumAmount2 +...
    LowerJejunumAmount2 + UpperIleumAmount2 + LowerIleumAmount2 ...
    ) / 1000;
    dglcLumenSIc = dglcLumenc2 - dglcLumenc1;

    cecumAmount1  = yout(i,308)*P_5737;
    colonAAmount1 = yout(i,317)*P_5772;
    colonTAmount1 = yout(i,326)*P_5807;
    colonDAmount1 = yout(i,335)*P_5840;
    colonSAmount1 = yout(i,344)*P_5877;
    dglcLumenLIc1 = (cecumAmount1 + colonAAmount1 + colonTAmount1 +...
        colonDAmount1 + colonSAmount1 ) / 1000;   
    cecumAmount2  = yout(i+1,308)*P_5737;
    colonAAmount2 = yout(i+1,317)*P_5772;
    colonTAmount2 = yout(i+1,326)*P_5807;
    colonDAmount2 = yout(i+1,335)*P_5840;
    colonSAmount2 = yout(i+1,344)*P_5877;
    dglcLumenLIc2 = (cecumAmount2 + colonAAmount2 + colonTAmount2 +...
        colonDAmount2 + colonSAmount2 ) / 1000;
    dglcLumenLIc = dglcLumenLIc2 - dglcLumenLIc1;
    %Stomach
    dglcLumenStomachc = (0-(IIf((EvalParameter(P_565, Time, y) & P_2743),...
        (EvalParameter(P_8156, Time, y)*y(258)*P_10749),0)))/1000;%y(258)

    %Feces
    dglcLumenFecesc = (IIf(EvalParameter(P_565, Time, y),(P_5911*y(353)),0))/1000;%y(362)
    %Add reactions
    %Transit stomach to SI 
    harvey = addReaction(harvey,'Stomach_SI_transit',{'glc_D[luSI]','glc_D[luSt]'},[1,-1],false,[],[],[],[],'',[],[],[]);
    %Stomach 
    harvey = addReaction(harvey,'Stomach_Lumen_uptake',{'Stomach_glc_D[c]','glc_D[luSt]'},[-1,1],true,[],[],[],[],'',[],[],[]);

    % fix b
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'glc_D[luLI]')),harvey.mets))) = dglcLumenSIc;% Lumne_glc_D[luLI]
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'glc_D[luSI]')),harvey.mets))) = dglcLumenLIc;% Lumne_glc_D[luSI]
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'glc_D[luSt]')),harvey.mets))) = dglcLumenStomachc*5;
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'glc_D[fe]')),harvey.mets)))   = dglcLumenFecesc*5;
    
    %print constraints
    dglcLumenSIc      = dglcLumenSIc
    dglcLumenLIc      = dglcLumenLIc
    dglcLumenStomachc = dglcLumenStomachc*5
    dglcLumenFecesc   = dglcLumenFecesc*5
    %%
    %look for Venous blood reactions
    %%
    %Constraints
    %Concentrations
    % Venous blood_glc_D(c)
    dglcbc = ((0-(P_4670*P_2754*((P_9835*P_5334)-(P_9833*(P_5341/P_7396))))) ...
            +(0-(((P_6016*(1-P_40))+P_1821)*P_5334)) ...
            +(P_75*P_5384) ...
            +(((P_3805*(1-P_40))-P_75)*P_5367) ...
            +(((P_3840*(1-P_40))-P_154)*P_5395) ...
            +(P_154*P_5413) ...
            +(P_230*P_5441) ...
            +(((P_3866*(1-P_40))-P_230)*P_5424) ...
            +(((P_3891*(1-P_40))-P_345)*P_5454) ...
            +(P_345*P_5471) ...
            +(((P_3918*(1-P_40))-P_418)*P_5482) ...
            +(P_418*P_5499) ...
            +(P_488*P_5529) ...
            +(((P_3945*(1-P_40))-P_488)*P_5511) ...
            +(P_825*P_5964) ...
            +(P_11064*P_11569) ...
            +(P_11232*P_11646) ...
            +(P_1673*P_6003) ...
            +((((P_4441+P_6109)*(1-P_40))-(P_1673+P_11708))*P_5986) ...
            +(P_1821*P_6035) ...
            +(P_1895*P_6063) ...
            +(((P_4493*(1-P_40))-P_1895)*P_6046) ...
            +(P_2001*P_6093) ...
            +(P_2118*P_6145) ...
            +(((P_4554*(1-P_40))-P_2118)*P_6128) ...
            +(P_2189*P_6173) ...
            +(P_10782*P_11152) ...
            +(P_10774*P_11120) ...
            +(P_10768*P_11088) ...
            +(P_10796*P_11216) ...
            +(P_10804*P_11255) ...
            +(P_10789*P_11184) ...
            +(P_10838*P_11405) ...
            +(P_10832*P_11375) ...
            +(P_10824*P_11345) ...
            +(P_10817*P_11315) ...
            +(P_10810*P_11285) ...
            +(((P_6185*(1-P_40))-P_6186)*P_8730) ...
            +(P_6186*P_8748) ...
            +(((P_4647*(1-P_40))-P_4648)*P_6196) ...
            +(P_4648*P_6210) ...
            +EvalParameter(P_3626, Time, y) ...
            +EvalParameter(P_3632, Time, y) ...
            +EvalParameter(P_3638, Time, y) ...
            +EvalParameter(P_3644, Time, y) ...
            +EvalParameter(P_3650, Time, y) ...
            +EvalParameter(P_3656, Time, y) ...
            +EvalParameter(P_3662, Time, y) ...
            +EvalParameter(P_3668, Time, y) ...
            +EvalParameter(P_3674, Time, y) ...
            +EvalParameter(P_3680, Time, y) ...
            +EvalParameter(P_3686, Time, y) ...
            +EvalParameter(P_3692, Time, y) ...
            +EvalParameter(P_3698, Time, y) ...
            +EvalParameter(P_3704, Time, y))/1000;%y(1)
    harvey.b(find(cellfun(@(x) ~isempty(strmatch(x,'glc_D[bc]')),harvey.mets))) = dglcbc*5;
    
    %print constraints
    dglcbc = dglcbc*5
    %%
%     Venous blood cells
%     findHarveyReactions('RBC',harvey);
%     findHarveyReactions('Platelet',harvey);
%     findHarveyReactions('Monocyte',harvey);
%     findHarveyReactions('Nkcells',harvey);
%     findHarveyReactions('CD4Tcells',harvey);
%     findHarveyReactions('Bcells',harvey);
    %%
    %Constraints (here all reactions that procude glc in cells are eq to dglc)
    %Concentrations
    % bloodcells_glc_D(c)
    dglcbloodcellc = ((0-P_9988) ...
            +(P_4670*P_2754*((P_9835*P_5334)-(P_9833*(P_5341/P_7396)))) ...
            +(0-(P_6016*P_40*P_5341)) ...
            +(P_3805*P_40*P_5375) ...
            +(P_3840*P_40*P_5404) ...
            +(P_3866*P_40*P_5432) ...
            +(P_3891*P_40*P_5462) ...
            +(P_3918*P_40*P_5490) ...
            +(P_3945*P_40*P_5520) ...
            +((P_4441+P_6109)*P_40*P_5994) ...
            +(P_4493*P_40*P_6054) ...
            +(P_4554*P_40*P_6136) ...
            +(P_6185*P_40*P_8737) ...
            +(P_4647*P_40*P_6203))/1000;%y(5)
    met='glc_D[c]';    
    [harvey]=metPool(harvey,met);
    harvey.b(end) = dglcbloodcellc*5;% bloodcells_g6p
    
    %bloodcells_g6p 
    met='g6p[c]';   
    dg6pbloodcellc = P_9988/1000;%y(7)
    [harvey]=metPool(harvey,met);
    harvey.b(end) = dg6pbloodcellc*5;% bloodcells_g6p

    %Hexokinase
    rxnName={'RBC_HEX1'...
        'Platelet_HEX1'...
        'Platelet_r0354'...
        'Platelet_r0355'...
        'Monocyte_HEX1'...
        'Monocyte_r0354'...
        'Monocyte_r0355'...
        'Nkcells_HEX1'...
        'Nkcells_r0354'...
        'Nkcells_r0355'...
        'CD4Tcells_HEX1'...
        'CD4Tcells_r0354'...
        'CD4Tcells_r0355'...
        'Bcells_HEX1'...
        'Bcells_r0354'...
        'Bcells_r0355'};
    oldNumVec = find(cellfun(@(x) ~isempty(strmatch(x,rxnName)),harvey.rxns));
    harvey = addMetabolite2Harvey(harvey,oldNumVec,[],[],ones(1,length(rxnName)));
    harvey.b(end) = P_9988/1000*5;
    
    %print constraints
    dglcbloodcellc = dglcbloodcellc*5
    dg6pbloodcellc = dg6pbloodcellc*5
    hexbloodcellc  = P_9988/1000*5