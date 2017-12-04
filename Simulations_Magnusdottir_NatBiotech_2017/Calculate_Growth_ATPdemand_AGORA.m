% This script calculates and plots the growth rates and the maximal ATP 
% production on the two diets used in Magnusdottir et al., Nat Biotech 2017.
% Growth rates are same as shown in Table S6, Magnusdottir et al., Nat Biotech 2017.
% Version of AGORA used: 1.0
% Stefania Magnusdottir & Almut Heinken, 01.12.2017

% initialize the COBRA Toolbox
initCobraToolbox

% Import file with information on AGORA including reconstruction names
currentDir = pwd;
[~,infoFile,~]=xlsread(strcat(currentDir,'\InputFiles\','ModelInformation.xlsx'));

%% do the calculations
Biomass.species=infoFile(2:end,3);
% Growth rates, Western diet
Biomass.WDanaerobic=zeros(size(Biomass.species));
Biomass.WDaerobic=zeros(size(Biomass.species));
% Growth rates, High fiber diet
Biomass.HFanaerobic=zeros(size(Biomass.species));
Biomass.HFaerobic=zeros(size(Biomass.species));

ATPDemand.species=infoFile(2:end,3);
% ATP demand, Western diet
ATPDemand.WDanaerobic=zeros(size(ATPDemand.species));
ATPDemand.WDaerobic=zeros(size(ATPDemand.species));
% ATP demand, High fiber diet
ATPDemand.HFanaerobic=zeros(size(ATPDemand.species));
ATPDemand.HFaerobic=zeros(size(ATPDemand.species));

% loop through all AGORA models and calculate the growth rates and ATP
% demand on the two diets with/without oxygen
for i=2:size(infoFile,1)
    load(strcat(currentDir,'\AGORA\',infoFile{i,3},'.mat'));
    % calculate growth rates
    bioID=model.rxns(strmatch('biomass',model.rxns));
    model=changeObjective(model,bioID{1,1});
    % Western diet
    model=useWesternDiet_AGORA(model);
    FBA=optimizeCbModel(model,'max');
    Biomass.WDanaerobic(i-1,1)=FBA.f;
    model=changeRxnBounds(model,'EX_o2(e)',-10,'l');
    FBA=optimizeCbModel(model,'max');
    Biomass.WDaerobic(i-1,1)=FBA.f;
    % High fiber diet
    model=useHighFiberDiet_AGORA(model);
    FBA=optimizeCbModel(model,'max');
    Biomass.HFanaerobic(i-1,1)=FBA.f;
    model=changeRxnBounds(model,'EX_o2(e)',-10,'l');
    FBA=optimizeCbModel(model,'max');
    Biomass.HFaerobic(i-1,1)=FBA.f;
    % calculate ATP demand
    model=changeObjective(model,'DM_atp_c_');
    % Western diet
    model=useWesternDiet_AGORA(model);
    FBA=optimizeCbModel(model,'max');
    ATPDemand.WDanaerobic(i-1,1)=FBA.f;
    model=changeRxnBounds(model,'EX_o2(e)',-10,'l');
    FBA=optimizeCbModel(model,'max');
    ATPDemand.WDaerobic(i-1,1)=FBA.f;
    % High fiber diet
    model=useHighFiberDiet_AGORA(model);
    FBA=optimizeCbModel(model,'max');
    ATPDemand.HFanaerobic(i-1,1)=FBA.f;
    model=changeRxnBounds(model,'EX_o2(e)',-10,'l');
    FBA=optimizeCbModel(model,'max');
    ATPDemand.HFaerobic(i-1,1)=FBA.f;
end

%% plot results
% plot the growth rates
dataAll=vertcat(Biomass.WDanaerobic,Biomass.WDaerobic,Biomass.HFanaerobic,Biomass.HFaerobic);
group={};
for i=1:length(Biomass.WDanaerobic)
    group{i,1}='Western diet, anoxic';
end
currSize=length(group);
for i=1:length(Biomass.WDaerobic)
    group{i+currSize,1}='Western diet, oxic';
end
currSize=length(group);
for i=1:length(Biomass.HFanaerobic)
    group{i+currSize,1}='High fiber diet, anoxic';
end
currSize=length(group);
for i=1:length(Biomass.HFaerobic)
    group{i+currSize,1}='High fiber diet, oxic';
end

figure;
boxplot(dataAll,group,'PlotStyle','traditional','BoxStyle','outline')
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    if j==1
        patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',.5);
    end
    if j==2
        patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',.5);
    end
    if j==3
        patch(get(h(j),'XData'),get(h(j),'YData'),'r','FaceAlpha',.5);
    end
    if j==4
        patch(get(h(j),'XData'),get(h(j),'YData'),'m','FaceAlpha',.5);
    end
end
set(gca,'FontSize',12,'XTickLabelRotation',90)
ylabel('Growth rate (hr-1)');
title('Growth rates calculated for AGORA models')
print('GrowthRates_AGORA','-dpng','-r300')

% plot the ATP demand
dataAll=vertcat(ATPDemand.WDanaerobic,ATPDemand.WDaerobic,ATPDemand.HFanaerobic,ATPDemand.HFaerobic);
group={};
for i=1:length(ATPDemand.WDanaerobic)
    group{i,1}='Western diet, anoxic';
end
currSize=length(group);
for i=1:length(ATPDemand.WDaerobic)
    group{i+currSize,1}='Western diet, oxic';
end
currSize=length(group);
for i=1:length(ATPDemand.HFanaerobic)
    group{i+currSize,1}='High fiber diet, anoxic';
end
currSize=length(group);
for i=1:length(ATPDemand.HFaerobic)
    group{i+currSize,1}='High fiber diet, oxic';
end

figure;
boxplot(dataAll,group,'PlotStyle','traditional','BoxStyle','outline')
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    if j==1
        patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',.5);
    end
    if j==2
        patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',.5);
    end
    if j==3
        patch(get(h(j),'XData'),get(h(j),'YData'),'r','FaceAlpha',.5);
    end
    if j==4
        patch(get(h(j),'XData'),get(h(j),'YData'),'m','FaceAlpha',.5);
    end
end
set(gca,'FontSize',12,'XTickLabelRotation',90)
ylabel('Flux through ATP demand reaction (mmol*gDW-1*hr-1)');
title('ATP production calculated for AGORA models')
print('ATPproduction_AGORA','-dpng','-r300')

%% create an output table

GrowthRates={'Model','WesternDiet_Anoxic','WesternDiet_Oxic','HighFiberDiet_Anoxic','FiberDiet_Oxic'};
for i=1:length(Biomass.species)
    GrowthRates{i+1,1}=Biomass.species{i};
    GrowthRates{i+1,2}=Biomass.WDanaerobic(i);
    GrowthRates{i+1,3}=Biomass.WDaerobic(i);
    GrowthRates{i+1,4}=Biomass.HFanaerobic(i);
    GrowthRates{i+1,5}=Biomass.HFaerobic(i);
end
save GrowthRates GrowthRates;

ATPproduction={'Model','WesternDiet_Anoxic','WesternDiet_Oxic','HighFiberDiet_Anoxic','FiberDiet_Oxic'};
for i=1:length(Biomass.species)
    ATPproduction{i+1,1}=ATPDemand.species{i};
    ATPproduction{i+1,2}=ATPDemand.WDanaerobic(i);
    ATPproduction{i+1,3}=ATPDemand.WDaerobic(i);
    ATPproduction{i+1,4}=ATPDemand.HFanaerobic(i);
    ATPproduction{i+1,5}=ATPDemand.HFaerobic(i);
end
save ATPproduction ATPproduction;


