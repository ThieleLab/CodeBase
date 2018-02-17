% Analysis of data performed in Heinken et al., Personalized modeling of the human 
% gut microbiome reveals distinct bile acid deconjugation and 
% biotransformation potential in healthy and IBD individuals (preprint on 
% bioRxiv, 2017)
% Almut Heinken, 01/2018

%% Execute this script to repeat the analysis.
% First, please download the AGORA resource (version 1.02) from here: https://webdav-r3lab.uni.lu/public/msp/AGORA-1.02/Agora-1.02-EuropeanAverage.zip
% Extract the mat files and the them into a folder.
% Define the path to the folder contaibning AGORA.
modelPath='YOUR_PATH_TO_AGORA';

% Import a file with information on the AGORA organisms including 
% reconstruction names and taxonomy.
[~,taxonomy,~]=xlsread('AGORA_infoFile.xlsx');

% Define the path where you want the created pairwise models and 
% computed results to  be saved.
savePath='YOUR_PATH_TO_RESULTS';

%% Calculation of reaction abundances

rxnsList={'TCDCHOLBHSe';'GCDCHOLBHSe';'GCHOLBHSe';'TCHOLBHSe';'12aHSDHe';'7AHSDHe';'CDCA7aHSDHe';'7AHSDH';'UCA7bHSDHe';'UDCA7bHSDHe';'CA3aHSDHe';'CDCA3aHSDHe';'ICA3bHSDHe';'ICDCA3bHSDHe';'BICoAL1';'BICoAL2';'BICoAL3';'BAIA1';'BAIA2';'BAIA3';'BAICDH1';'BAICDH2';'BAICDH3';'BAICDH4';'BAICDH5';'BAICDH6';'BAIF1';'BAIF2';'BAIF3';'BAIEI1';'BAIEI2';'BAIEI3'};

% calculation of reaction abundance for HMP microbiomes
load('HMP_abundance.mat');
[ReactionAbundance]=calculateReactionAbundance(HMP_abundance,modelPath,taxonomy,rxnsList,[]);
save(strcat(savePath,'HMP_ReactionAbundance'),'ReactionAbundance');
%delete(gcp('nocreate'))
% calculation of reaction abundance for MetaHIT microbiomes
load('MetaHIT_abundance.mat');
[ReactionAbundance]=calculateReactionAbundance(MetaHIT_abundance,modelPath,taxonomy,rxnsList,[]);
save(strcat(savePath,'MetaHIT_ReactionAbundance'),'ReactionAbundance');
%delete(gcp('nocreate'))
% calculation of reaction abundance for pediatric IBD microbiomes
load('pIBD_abundance.mat');
[ReactionAbundance]=calculateReactionAbundance(pIBD_abundance,modelPath,taxonomy,rxnsList,[]);
save(strcat(savePath,'pIBD_ReactionAbundance'),'ReactionAbundance');


%% Creation of output tables
%% Prediction of bile acid production by single models
% Initialize the COBRA Toolbox
initCobraToolbox
% create a table with bile acid production by the 217 bile acid producers
% on Average European diet
load('BileAcidProducers.mat');
BA_Objectives={'EX_C02528(e)','EX_cholate(e)','EX_12dhchol(e)','EX_7ocholate(e)','EX_7dhcdchol(e)','EX_3dhchol(e)','EX_3dhcdchol(e)','EX_isochol(e)','EX_icdchol(e)','EX_HC02191(e)','EX_dchac(e)','EX_adchac(e)','EX_alchac(e)','EX_uchol(e)','EX_HC02194(e)'};
for j=1:length(BA_Objectives)
    SingleModels_Output_Table{1,j+1}=BA_Objectives{j};
end
for i=1:length(BileAcidProducers)
    load(strcat(modelPath,BileAcidProducers{i,1}));
    SingleModels_Output_Table{i+1,1}=BileAcidProducers{i,1};
    for j=1:length(BA_Objectives)
        if ~isempty(find(strcmp(model.rxns,BA_Objectives{j})))
            SingleModels_Output_Table{1,j+1}=model.rxnNames(find(strcmp(model.rxns,BA_Objectives{j})));
            model=changeObjective(model,BA_Objectives{j});
            FBAsolution=optimizeCbModel(model,'max');
            SingleModels_Output_Table{i+1,j+1}=FBAsolution.f;
        end
    end
end
save(strcat(savePath,'SingleModels_Output_Table'),'SingleModels_Output_Table');

%% Analysis of pairwise models
load('BileAcid_Production_ComplementaryPairs.mat');
% calculate the number of pairs that enabled production for each bile acid
mets=unique(ComplementaryPairs(2:end,9));
for i=1:length(mets)
    Number_Enabling_Pairs{i,1}=mets{i};
    numPairs=strmatch(mets{i},ComplementaryPairs(:,9));
    Number_Enabling_Pairs{i,2}=size(numPairs,1);
end

% calculate the number of total bile acids produced by each single model and pair (shown in
% Figure 3a)
load('Pairwise_BileAcid_Production.mat');
ProducerPairs={};
for i=1:length(BileAcidProducers)
    ProducerPairs{i+1,1}=BileAcidProducers{i,1};
    for j=1:length(BileAcidProducers)
        ProducerPairs{1,j+1}=BileAcidProducers{j,1};
        ProducerPairs{i+1,j+1}=0;
        str1InRow1=find(strcmp(Pairwise_BileAcid_Production(:,1),BileAcidProducers{i,1}));
        str2InRow2=find(strcmp(Pairwise_BileAcid_Production(:,2),BileAcidProducers{j,1}));
        % find the microbe pairs in both rows
        findpairs=intersect(str1InRow1,str2InRow2);
        if ~isempty(findpairs)
            % length=number of cases of complementary production
            ProducerPairs{i+1,j+1}=ProducerPairs{i+1,j+1}+length(findpairs);
        end
        if ~strcmp(BileAcidProducers{i,1},BileAcidProducers{j,1})
            str1InRow2=find(strcmp(Pairwise_BileAcid_Production(:,2),BileAcidProducers{i,1}));
            str2InRow1=find(strcmp(Pairwise_BileAcid_Production(:,1),BileAcidProducers{j,1}));
            % find the microbe pairs in both rows
            findpairs=intersect(str1InRow2,str2InRow1);
            if ~isempty(findpairs)
                % length=number of cases of complementary production
                ProducerPairs{i+1,j+1}=ProducerPairs{i+1,j+1}+length(findpairs);
            end
        end
    end
end
save(strcat(savePath,'BileAcid_Production_ProducerPairs_forFigure3a'),'ProducerPairs');

%% analysis of HMP microbiome models
load('HMP_IndividualIDs.mat');
load('HMP_BileAcid_Production.mat');
% summarize the computed bile acid production potential values (data shown
% in Table S6)
mets=fieldnames(HMP_BileAcid_Production);
mets=sort(mets);
for i=1:length(HMP_IndividualIDs)
    HMP_Output_Table{i+1,1}=HMP_IndividualIDs{i,1};
    for j=1:length(mets)
        HMP_Output_Table{1,j+1}=strrep(mets{j},'BA_','');
        if ~isempty(HMP_BileAcid_Production.(mets{j})(strcmp(HMP_IndividualIDs{i,1},HMP_BileAcid_Production.(mets{j})(:,1)),2))
            val=HMP_BileAcid_Production.(mets{j})(strcmp(HMP_IndividualIDs{i,1},HMP_BileAcid_Production.(mets{j})(:,1)),2);
            val=val{1};
            if ischar(val)
                HMP_Output_Table{i+1,j+1}=str2num(val);
            else
                HMP_Output_Table{i+1,j+1}=val;
            end
        else
            HMP_Output_Table{i+1,j+1}=0;
        end
    end
end
save(strcat(savePath,'HMP_Output_Table'),'HMP_Output_Table');
% normalize the computed bile acid production potential values (data shown
% in Figure 3b)
for i=1:length(HMP_IndividualIDs)
    HMP_Output_Table_Normalized{i+1,1}=HMP_IndividualIDs{i,1};
    for j=1:length(mets)
        % determine the highest flux that was achieved
        for k=2:length(HMP_BileAcid_Production.(mets{j}))
            if ~isempty(HMP_BileAcid_Production.(mets{j}){k,2})
                allVal(k-1,1)=str2num(HMP_BileAcid_Production.(mets{j}){k,2});
            else
                allVal(k-1,1)=0;
            end
        end
        HMP_Output_Table_Normalized{1,j+1}=strrep(mets{j},'BA_','');
        topval=max(allVal);
        if topval==0
            HMP_Output_Table_Normalized{i+1,j+1}=0;
        else
            if ~isempty(HMP_BileAcid_Production.(mets{j})(strcmp(HMP_IndividualIDs{i,1},HMP_BileAcid_Production.(mets{j})(:,1)),2))
                val=HMP_BileAcid_Production.(mets{j})(strcmp(HMP_IndividualIDs{i,1},HMP_BileAcid_Production.(mets{j})(:,1)),2);
                val=val{1};
                if ischar(val)
                    val=str2num(val);
                end
                HMP_Output_Table_Normalized{i+1,j+1}=val\topval;
            else
                HMP_Output_Table_Normalized{i+1,j+1}=0;
            end
        end
    end
end
save(strcat(savePath,'HMP_Output_Table_Normalized_forFigure3b'),'HMP_Output_Table_Normalized');

% load the flux spans on internal exchanges representing strain-level
% contributions to overall production (data shown in Figure S4)
load('HMP_FluxSpans_InternalExchanges.mat');

%% analysis of MetaHIT microbiome models
load('MetaHIT_IndividualIDs.mat');
load('MetaHIT_BileAcid_Production.mat');
% summarize the computed bile acid production potential values (data shown
% in Table S7)
mets=fieldnames(MetaHIT_BileAcid_Production);
mets=sort(mets);
for i=1:length(MetaHIT_IndividualIDs)
    MetaHIT_Output_Table{i+1,1}=MetaHIT_IndividualIDs{i,1};
    for j=1:length(mets)
        MetaHIT_Output_Table{1,j+1}=strrep(mets{j},'BA_','');
        if ~isempty(MetaHIT_BileAcid_Production.(mets{j})(strcmp(MetaHIT_IndividualIDs{i,1},MetaHIT_BileAcid_Production.(mets{j})(:,1)),2))
            val=MetaHIT_BileAcid_Production.(mets{j})(strcmp(MetaHIT_IndividualIDs{i,1},MetaHIT_BileAcid_Production.(mets{j})(:,1)),2);
            val=val{1};
            if ischar(val)
                MetaHIT_Output_Table{i+1,j+1}=str2num(val);
            else
                MetaHIT_Output_Table{i+1,j+1}=val;
            end
        else
            MetaHIT_Output_Table{i+1,j+1}=0;
        end
    end
end
save(strcat(savePath,'MetaHIT_Output_Table'),'MetaHIT_Output_Table');
% normalize the computed bile acid production potential values (data shown
% in Figure 3c)
for i=1:length(MetaHIT_IndividualIDs)
    MetaHIT_Output_Table_Normalized{i+1,1}=MetaHIT_IndividualIDs{i,1};
    for j=1:length(mets)
        % determine the highest flux that was achieved
        for k=2:length(MetaHIT_BileAcid_Production.(mets{j}))
            if ~isempty(MetaHIT_BileAcid_Production.(mets{j}){k,2})
                allVal(k-1,1)=str2num(MetaHIT_BileAcid_Production.(mets{j}){k,2});
            else
                allVal(k-1,1)=0;
            end
        end
        MetaHIT_Output_Table_Normalized{1,j+1}=strrep(mets{j},'BA_','');
        topval=max(allVal);
        if topval==0
            MetaHIT_Output_Table_Normalized{i+1,j+1}=0;
        else
            if ~isempty(MetaHIT_BileAcid_Production.(mets{j})(strcmp(MetaHIT_IndividualIDs{i,1},MetaHIT_BileAcid_Production.(mets{j})(:,1)),2))
                val=MetaHIT_BileAcid_Production.(mets{j})(strcmp(MetaHIT_IndividualIDs{i,1},MetaHIT_BileAcid_Production.(mets{j})(:,1)),2);
                val=val{1};
                if ischar(val)
                    val=str2num(val);
                end
                MetaHIT_Output_Table_Normalized{i+1,j+1}=val\topval;
            else
                MetaHIT_Output_Table_Normalized{i+1,j+1}=0;
            end
        end
    end
end
save(strcat(savePath,'MetaHIT_Output_Table_Normalized_forFigure3c'),'MetaHIT_Output_Table_Normalized');

% load the flux spans on internal exchanges representing strain-level
% contributions to overall production (data shown in Figure S5)
load('MetaHIT_FluxSpans_InternalExchanges.mat');

% analysis of pediatric IBD microbiome models
load('pIBD_IndividualIDs.mat');
load('pIBD_BileAcid_Production.mat');
% summarize the computed bile acid production potential values (data shown
% in Table S8)
mets=fieldnames(pIBD_BileAcid_Production);
mets=sort(mets);
for i=1:length(pIBD_IndividualIDs)
    pIBD_Output_Table{i+1,1}=pIBD_IndividualIDs{i,1};
    for j=1:length(mets)
        pIBD_Output_Table{1,j+1}=strrep(mets{j},'BA_','');
        if ~isempty(pIBD_BileAcid_Production.(mets{j})(strcmp(pIBD_IndividualIDs{i,1},pIBD_BileAcid_Production.(mets{j})(:,1)),2))
            val=pIBD_BileAcid_Production.(mets{j})(strcmp(pIBD_IndividualIDs{i,1},pIBD_BileAcid_Production.(mets{j})(:,1)),2);
            val=val{1};
            if ischar(val)
                pIBD_Output_Table{i+1,j+1}=str2num(val);
            else
                pIBD_Output_Table{i+1,j+1}=val;
            end
        else
            pIBD_Output_Table{i+1,j+1}=0;
        end
    end
end
save(strcat(savePath,'pIBD_Output_Table'),'pIBD_Output_Table');
% normalize the computed bile acid production potential values (data shown
% in Figure 3d)
for i=1:length(pIBD_IndividualIDs)
    pIBD_Output_Table_Normalized{i+1,1}=pIBD_IndividualIDs{i,1};
    for j=1:length(mets)
        % determine the highest flux that was achieved
        for k=2:length(pIBD_BileAcid_Production.(mets{j}))
            if ~isempty(pIBD_BileAcid_Production.(mets{j}){k,2})
                allVal(k-1,1)=str2num(pIBD_BileAcid_Production.(mets{j}){k,2});
            else
                allVal(k-1,1)=0;
            end
        end
        pIBD_Output_Table_Normalized{1,j+1}=strrep(mets{j},'BA_','');
        topval=max(allVal);
        if topval==0
            pIBD_Output_Table_Normalized{i+1,j+1}=0;
        else
            if ~isempty(pIBD_BileAcid_Production.(mets{j})(strcmp(pIBD_IndividualIDs{i,1},pIBD_BileAcid_Production.(mets{j})(:,1)),2))
                val=pIBD_BileAcid_Production.(mets{j})(strcmp(pIBD_IndividualIDs{i,1},pIBD_BileAcid_Production.(mets{j})(:,1)),2);
                val=val{1};
                if ischar(val)
                    val=str2num(val);
                end
                pIBD_Output_Table_Normalized{i+1,j+1}=val\topval;
            else
                pIBD_Output_Table_Normalized{i+1,j+1}=0;
            end
        end
    end
end
save(strcat(savePath,'pIBD_Output_Table_Normalized_forFigure3d'),'pIBD_Output_Table_Normalized');

% load the flux spans on internal exchanges representing strain-level
% contributions to overall production (data shown in Figure 4a)
load('pIBD_FluxSpans_InternalExchanges.mat');

%% Statistical analysis
% MetaHIT models: comparison of healthy and ulcerative colitis microbiomes
% pIBD models: comparison of healthy and pediatric Crohn's Disease
% microbiomes

%% MetaHIT dataset
% bile acid production quantitative
load('MetaHIT_IndividualIDs.mat');
load(strcat(savePath,'MetaHIT_Output_Table.mat'));
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_UC';
Statistics{1,8}='SD_UC';

for i=2:size(MetaHIT_Output_Table,2)
    Statistics{i,1}=MetaHIT_Output_Table{1,i};
    data1=[];
    data2=[];
    cnt=1;
    for j=2:size(MetaHIT_Output_Table,1)
        if strcmp(MetaHIT_IndividualIDs{find(strcmp(MetaHIT_Output_Table{j,1},MetaHIT_IndividualIDs(:,1))),2},'healthy')
            data1(cnt,1)=MetaHIT_Output_Table{j,i};
            cnt=cnt+1;
        end
    end
    cnt=1;
    for j=2:size(MetaHIT_Output_Table,1)
        if strcmp(MetaHIT_IndividualIDs{find(strcmp(MetaHIT_Output_Table{j,1},MetaHIT_IndividualIDs(:,1))),2},'UC')
            data2(cnt,1)=MetaHIT_Output_Table{j,i};
            cnt=cnt+1;
        end
    end
    [p,h,stats] = ranksum(data1,data2);
    Statistics{i,2}=p;
    Statistics{i,3}=h;
    Statistics{i,4}=stats;
    Statistics{i,5}=mean(data1);
    Statistics{i,6}=std(data1);
    Statistics{i,7}=mean(data2);
    Statistics{i,8}=std(data2);
end
save(strcat(savePath,'MetaHIT_Statistics_BA_Production'),'Statistics');

% total reaction abundance
load(strcat(savePath,'MetaHIT_ReactionAbundance.mat'));
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_UC';
Statistics{1,8}='SD_UC';

for i=2:size(ReactionAbundance.('Total'),2)
    Statistics{i,1}=ReactionAbundance.('Total'){1,i};
    data1=[];
    data2=[];
    if ~isempty(ReactionAbundance.('Total'){i,2})
        cnt=1;
        for j=2:size(ReactionAbundance.('Total'),1)
            if strcmp(MetaHIT_IndividualIDs{find(strcmp(ReactionAbundance.('Total'){j,1},MetaHIT_IndividualIDs(:,1))),2},'healthy')
                data1(cnt,1)=ReactionAbundance.('Total'){j,i};
                cnt=cnt+1;
            end
        end
        cnt=1;
        for j=2:size(ReactionAbundance.('Total'),1)
            if strcmp(MetaHIT_IndividualIDs{find(strcmp(ReactionAbundance.('Total'){j,1},MetaHIT_IndividualIDs(:,1))),2},'UC')
                data2(cnt,1)=ReactionAbundance.('Total'){j,i};
                cnt=cnt+1;
            end
        end
        [p,h,stats] = ranksum(data1,data2);
        Statistics{i,2}=p;
        Statistics{i,3}=h;
        Statistics{i,4}=stats;
        Statistics{i,5}=mean(data1);
        Statistics{i,6}=std(data1);
        Statistics{i,7}=mean(data2);
        Statistics{i,8}=std(data2);
    end
end
save(strcat(savePath,'MetaHIT_Statistics_TotalReactionAbundance'),'Statistics');

% Reaction abundance on the phylum level
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_UC';
Statistics{1,8}='SD_UC';

for i=2:size(ReactionAbundance.('Phylum'),2)
    Statistics{i,1}=ReactionAbundance.('Phylum'){1,i};
    data1=[];
    data2=[];
    if ~isempty(ReactionAbundance.('Phylum'){2,i})
        cnt=1;
        for j=2:size(ReactionAbundance.('Phylum'),1)
            if strcmp(MetaHIT_IndividualIDs{find(strcmp(ReactionAbundance.('Phylum'){j,1},MetaHIT_IndividualIDs(:,1))),2},'healthy')
                data1(cnt,1)=ReactionAbundance.('Phylum'){j,i};
                cnt=cnt+1;
            end
        end
        cnt=1;
        for j=2:size(ReactionAbundance.('Phylum'),1)
            if strcmp(MetaHIT_IndividualIDs{find(strcmp(ReactionAbundance.('Phylum'){j,1},MetaHIT_IndividualIDs(:,1))),2},'UC')
                data2(cnt,1)=ReactionAbundance.('Phylum'){j,i};
                cnt=cnt+1;
            end
        end
        [p,h,stats] = ranksum(data1,data2);
        Statistics{i,2}=p;
        Statistics{i,3}=h;
        Statistics{i,4}=stats;
        Statistics{i,5}=mean(data1);
        Statistics{i,6}=std(data1);
        Statistics{i,7}=mean(data2);
        Statistics{i,8}=std(data2);
    end
end
save(strcat(savePath,'MetaHIT_Statistics_ReactionAbundance_Phylum'),'Statistics');

% Reaction abundance on the genus level
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_UC';
Statistics{1,8}='SD_UC';

for i=2:size(ReactionAbundance.('Genus'),2)
    Statistics{i,1}=ReactionAbundance.('Genus'){1,i};
    data1=[];
    data2=[];
    if ~isempty(ReactionAbundance.('Genus'){2,i})
        cnt=1;
        for j=2:size(ReactionAbundance.('Genus'),1)
            if strcmp(MetaHIT_IndividualIDs{find(strcmp(ReactionAbundance.('Genus'){j,1},MetaHIT_IndividualIDs(:,1))),2},'healthy')
                data1(cnt,1)=ReactionAbundance.('Genus'){j,i};
                cnt=cnt+1;
            end
        end
        cnt=1;
        for j=2:size(ReactionAbundance.('Genus'),1)
            if strcmp(MetaHIT_IndividualIDs{find(strcmp(ReactionAbundance.('Genus'){j,1},MetaHIT_IndividualIDs(:,1))),2},'UC')
                data2(cnt,1)=ReactionAbundance.('Genus'){j,i};
                cnt=cnt+1;
            end
        end
        [p,h,stats] = ranksum(data1,data2);
        Statistics{i,2}=p;
        Statistics{i,3}=h;
        Statistics{i,4}=stats;
        Statistics{i,5}=mean(data1);
        Statistics{i,6}=std(data1);
        Statistics{i,7}=mean(data2);
        Statistics{i,8}=std(data2);
    end
end
save(strcat(savePath,'MetaHIT_Statistics_ReactionAbundance_Genus'),'Statistics');

% flux spans
load('MetaHIT_FluxSpans_InternalExchanges.mat');
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_UC';
Statistics{1,8}='SD_UC';

for i=4:size(MetaHIT_FluxSpans_InternalExchanges,1)
    Statistics{i,1}=MetaHIT_FluxSpans_InternalExchanges{i,1};
    data1=[];
    data2=[];
    if ~isempty(MetaHIT_FluxSpans_InternalExchanges{i,2})
        cnt=1;
        for j=4:size(MetaHIT_FluxSpans_InternalExchanges,2)
            if strcmp(MetaHIT_IndividualIDs{strcmp(MetaHIT_FluxSpans_InternalExchanges{1,j},MetaHIT_IndividualIDs(:,1)),2},'healthy')
                data1(cnt,1)=MetaHIT_FluxSpans_InternalExchanges{i,j};
                cnt=cnt+1;
            end
        end
        cnt=1;
        for j=4:size(MetaHIT_FluxSpans_InternalExchanges,2)
            if strcmp(MetaHIT_IndividualIDs{find(strcmp(MetaHIT_FluxSpans_InternalExchanges{1,j},MetaHIT_IndividualIDs(:,1))),2},'UC')
                data2(cnt,1)=MetaHIT_FluxSpans_InternalExchanges{i,j};
                cnt=cnt+1;
            end
        end
        [p,h,stats] = ranksum(data1,data2);
        Statistics{i,2}=p;
        Statistics{i,3}=h;
        Statistics{i,4}=stats;
        Statistics{i,5}=mean(data1);
        Statistics{i,6}=std(data1);
        Statistics{i,7}=mean(data2);
        Statistics{i,8}=std(data2);
    end
end
save(strcat(savePath,'MetaHIT_Statistics_FluxSpans'),'Statistics');

%% pediatric IBD dataset
% bile acid production quantitative
load('pIBD_IndividualIDs.mat');
load(strcat(savePath,'pIBD_Output_Table.mat'));
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_CD';
Statistics{1,8}='SD_CD';

for i=2:size(pIBD_Output_Table,2)
    Statistics{i,1}=pIBD_Output_Table{1,i};
    data1=[];
    data2=[];
    cnt=1;
    for j=2:size(pIBD_Output_Table,1)
        if strcmp(pIBD_IndividualIDs{find(strcmp(pIBD_Output_Table{j,1},pIBD_IndividualIDs(:,1))),2},'healthy')
            data1(cnt,1)=pIBD_Output_Table{j,i};
            cnt=cnt+1;
        end
    end
    cnt=1;
    for j=2:size(pIBD_Output_Table,1)
        if strcmp(pIBD_IndividualIDs{find(strcmp(pIBD_Output_Table{j,1},pIBD_IndividualIDs(:,1))),2},'CD')
            data2(cnt,1)=pIBD_Output_Table{j,i};
            cnt=cnt+1;
        end
    end
    [p,h,stats] = ranksum(data1,data2);
    Statistics{i,2}=p;
    Statistics{i,3}=h;
    Statistics{i,4}=stats;
    Statistics{i,5}=mean(data1);
    Statistics{i,6}=std(data1);
    Statistics{i,7}=mean(data2);
    Statistics{i,8}=std(data2);
end
save(strcat(savePath,'pIBD_Statistics_BA_Production'),'Statistics');

% total reaction abundance
load(strcat(savePath,'pIBD_ReactionAbundance.mat'));
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_CD';
Statistics{1,8}='SD_CD';

for i=2:size(ReactionAbundance.('Total'),2)
    Statistics{i,1}=ReactionAbundance.('Total'){1,i};
    data1=[];
    data2=[];
    if ~isempty(ReactionAbundance.('Total'){2,i})
        cnt=1;
        for j=2:size(ReactionAbundance.('Total'),1)
            if strcmp(pIBD_IndividualIDs{find(strcmp(ReactionAbundance.('Total'){j,1},pIBD_IndividualIDs(:,1))),2},'healthy')
                data1(cnt,1)=ReactionAbundance.('Total'){j,i};
                cnt=cnt+1;
            end
        end
        cnt=1;
        for j=2:size(ReactionAbundance.('Total'),1)
            if strcmp(pIBD_IndividualIDs{find(strcmp(ReactionAbundance.('Total'){j,1},pIBD_IndividualIDs(:,1))),2},'CD')
                data2(cnt,1)=ReactionAbundance.('Total'){j,i};
                cnt=cnt+1;
            end
        end
        [p,h,stats] = ranksum(data1,data2);
        Statistics{i,2}=p;
        Statistics{i,3}=h;
        Statistics{i,4}=stats;
        Statistics{i,5}=mean(data1);
        Statistics{i,6}=std(data1);
        Statistics{i,7}=mean(data2);
        Statistics{i,8}=std(data2);
    end
end
save(strcat(savePath,'pIBD_Statistics_TotalReactionAbundance'),'Statistics');

% Reaction abundance on the phylum level
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_CD';
Statistics{1,8}='SD_CD';

for i=2:size(ReactionAbundance.('Phylum'),2)
    Statistics{i,1}=ReactionAbundance.('Phylum'){1,i};
    data1=[];
    data2=[];
    if ~isempty(ReactionAbundance.('Phylum'){2,i})
        cnt=1;
        for j=2:size(ReactionAbundance.('Phylum'),1)
            if strcmp(pIBD_IndividualIDs{find(strcmp(ReactionAbundance.('Phylum'){j,1},pIBD_IndividualIDs(:,1))),2},'healthy')
                data1(cnt,1)=ReactionAbundance.('Phylum'){j,i};
                cnt=cnt+1;
            end
        end
        cnt=1;
        for j=2:size(ReactionAbundance.('Phylum'),1)
            if strcmp(pIBD_IndividualIDs{find(strcmp(ReactionAbundance.('Phylum'){j,1},pIBD_IndividualIDs(:,1))),2},'CD')
                data2(cnt,1)=ReactionAbundance.('Phylum'){j,i};
                cnt=cnt+1;
            end
        end
        [p,h,stats] = ranksum(data1,data2);
        Statistics{i,2}=p;
        Statistics{i,3}=h;
        Statistics{i,4}=stats;
        Statistics{i,5}=mean(data1);
        Statistics{i,6}=std(data1);
        Statistics{i,7}=mean(data2);
        Statistics{i,8}=std(data2);
    end
end
save(strcat(savePath,'pIBD_Statistics_ReactionAbundance_Phylum'),'Statistics');

% Reaction abundance on the genus level
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_CD';
Statistics{1,8}='SD_CD';

for i=2:size(ReactionAbundance.('Genus'),2)
    Statistics{i,1}=ReactionAbundance.('Genus'){1,i};
    data1=[];
    data2=[];
    if ~isempty(ReactionAbundance.('Genus'){2,i})
        cnt=1;
        for j=2:size(ReactionAbundance.('Genus'),1)
            if strcmp(pIBD_IndividualIDs{find(strcmp(ReactionAbundance.('Genus'){j,1},pIBD_IndividualIDs(:,1))),2},'healthy')
                data1(cnt,1)=ReactionAbundance.('Genus'){j,i};
                cnt=cnt+1;
            end
        end
        cnt=1;
        for j=2:size(ReactionAbundance.('Genus'),1)
            if strcmp(pIBD_IndividualIDs{find(strcmp(ReactionAbundance.('Genus'){j,1},pIBD_IndividualIDs(:,1))),2},'CD')
                data2(cnt,1)=ReactionAbundance.('Genus'){j,i};
                cnt=cnt+1;
            end
        end
        [p,h,stats] = ranksum(data1,data2);
        Statistics{i,2}=p;
        Statistics{i,3}=h;
        Statistics{i,4}=stats;
        Statistics{i,5}=mean(data1);
        Statistics{i,6}=std(data1);
        Statistics{i,7}=mean(data2);
        Statistics{i,8}=std(data2);
    end
end
save(strcat(savePath,'pIBD_Statistics_ReactionAbundance_Genus'),'Statistics');

% flux spans
load('pIBD_FluxSpans_InternalExchanges.mat');
Statistics={};
Statistics{1,1}='Reaction';
Statistics{1,2}='p_value';
Statistics{1,3}='Decision';
Statistics{1,4}='Test_statistic';
Statistics{1,5}='Average_Healthy';
Statistics{1,6}='SD_Healthy';
Statistics{1,7}='Average_CD';
Statistics{1,8}='SD_CD';

for i=4:size(pIBD_FluxSpans_InternalExchanges,1)
    Statistics{i,1}=pIBD_FluxSpans_InternalExchanges{i,1};
    data1=[];
    data2=[];
    if ~isempty(pIBD_FluxSpans_InternalExchanges{i,2})
        cnt=1;
        for j=4:size(pIBD_FluxSpans_InternalExchanges,2)
            if strcmp(pIBD_IndividualIDs{strcmp(pIBD_FluxSpans_InternalExchanges{1,j},pIBD_IndividualIDs(:,1)),2},'healthy')
                data1(cnt,1)=pIBD_FluxSpans_InternalExchanges{i,j};
                cnt=cnt+1;
            end
        end
        cnt=1;
        for j=4:size(pIBD_FluxSpans_InternalExchanges,2)
            if strcmp(pIBD_IndividualIDs{find(strcmp(pIBD_FluxSpans_InternalExchanges{1,j},pIBD_IndividualIDs(:,1))),2},'CD')
                data2(cnt,1)=pIBD_FluxSpans_InternalExchanges{i,j};
                cnt=cnt+1;
            end
        end
        [p,h,stats] = ranksum(data1,data2);
        Statistics{i,2}=p;
        Statistics{i,3}=h;
        Statistics{i,4}=stats;
        Statistics{i,5}=mean(data1);
        Statistics{i,6}=std(data1);
        Statistics{i,7}=mean(data2);
        Statistics{i,8}=std(data2);
    end
end
save(strcat(savePath,'pIBD_Statistics_FluxSpans'),'Statistics');

%% Analysis of shadow prices
% load shadow prices for pediatric IBD models to calculate the statistics
load('pIBD_ShadowPricesExtracted.mat');
% extract and summarize in table-data shown in Figure 5a
pIBD_ShadowPrices_Table={};
for i=2:size(ShadowPricesExtracted,1)
    pIBD_ShadowPrices_Table{1,i}=ShadowPricesExtracted{i,1};
end
rowCnt=2;
for j=2:size(ShadowPricesExtracted,2)
    SPmets={};
    cnt=1;
    for i=2:size(ShadowPricesExtracted,1)
        if ~isempty(ShadowPricesExtracted{i,j})
        for k=1:size(ShadowPricesExtracted{i,j},1)
            SPmets{cnt,1}=ShadowPricesExtracted{i,j}{k,1};
            cnt=cnt+1;
        end
        end
    end
    SPmets=unique(SPmets);
    for k=1:length(SPmets)
        pIBD_ShadowPrices_Table{rowCnt,1}=strcat(SPmets{k},'_',ShadowPricesExtracted{1,j});
        for i=2:size(ShadowPricesExtracted,1)
            if ~isempty(ShadowPricesExtracted{i,j})
                findMet=find(strcmp(ShadowPricesExtracted{i,j}(:,1),SPmets{k}));
                if ~isempty(findMet)
                    pIBD_ShadowPrices_Table{rowCnt,i}=ShadowPricesExtracted{i,j}{findMet,2};
                else
                    pIBD_ShadowPrices_Table{rowCnt,i}=0;
                end
            else
                pIBD_ShadowPrices_Table{rowCnt,i}=0;
            end
        end
        rowCnt=rowCnt+1;
    end
end
save(strcat(savePath,'pIBD_ShadowPrices_Table'),'pIBD_ShadowPrices_Table');

% calculate the statistics, comparing pediatric IBD patients in the dataset
% and healthy controls
% count the negative entries
for i=2:size(pIBD_ShadowPrices_Table,2)
    SP_Sums{i-1,1}=pIBD_ShadowPrices_Table{1,i};
    SP_Sums{i-1,2}=0;
    for j=2:size(pIBD_ShadowPrices_Table,1)
        if pIBD_ShadowPrices_Table{j,i}<0
        SP_Sums{i-1,2}=SP_Sums{i-1,2}+1;
        end
    end
end

cntH=1;
cntCD=1;
for i=1:length(SP_Sums)
    group=pIBD_IndividualIDs(find(strcmp(SP_Sums{i,1},pIBD_IndividualIDs(:,1))),2);
    if strcmp(group{1},'healthy')
        dataH(cntH,1)=SP_Sums{i,2};
        cntH=cntH+1;
    end
    if strcmp(group{1},'CD')
        dataCD(cntCD,1)=SP_Sums{i,2};
        cntCD=cntCD+1;
    end
end
[p,h,stats] = ranksum(dataH,dataCD);
Statistics{1,1}='p value';
Statistics{2,1}='Decision';
Statistics{3,1}='z-value';
Statistics{4,1}='Rank sum';
Statistics{5,1}='Mean healthy';
Statistics{6,1}='SD healthy';
Statistics{7,1}='Mean CD';
Statistics{8,1}='SD CD';
Statistics{1,2}=p;
Statistics{2,2}=h;
Statistics{3,2}=stats.zval;
Statistics{4,2}=stats.ranksum;
Statistics{5,2}=mean(dataH);
Statistics{6,2}=std(dataH);
Statistics{7,2}=mean(dataCD);
Statistics{8,2}=std(dataCD);

%% analysis of shadow prices in all microbiome models

% correlation of reaction abundances with bile acid production potential in
% all microbiome models (data shown in Table S11)
% reaction abundances are calculated from the mapped relative strain-level abundances
% and the reaction content of the corresponding AGORA models (see above)
load('pIBD_BileAcid_Production.mat');
mets=fieldnames(pIBD_BileAcid_Production);
for i=1:length(mets)
    %% first combine all three datasets
    load('pIBD_IndividualIDs.mat');
    load('pIBD_BileAcid_Production.mat');
    Production={};
    % production
    Production{1,1}='Individual';
    Production{1,2}='Stratification';
    Production{1,3}='Flux';
    cnt=2;
    for j=2:size(pIBD_BileAcid_Production.(mets{i,1}),1)
        Production{cnt,1}=pIBD_BileAcid_Production.(mets{i,1}){j,1};
        Production{cnt,3}=str2num(pIBD_BileAcid_Production.(mets{i,1}){j,2});
        group=pIBD_IndividualIDs(find(strcmp(pIBD_BileAcid_Production.(mets{i,1}){j,1},pIBD_IndividualIDs(:,1))),2);
        if strcmp(group{1},'healthy')
            Production{cnt,2}='Healthy_pIBD';
        end
        if strcmp(group{1},'CD')
            Production{cnt,2}='CD_pIBD';
        end
        cnt=cnt+1;
    end
    load('MetaHIT_BileAcid_Production.mat');
    load('MetaHIT_IndividualIDs.mat');
    for j=2:size(MetaHIT_BileAcid_Production.(mets{i,1}),1)
        Production{cnt,1}=MetaHIT_BileAcid_Production.(mets{i,1}){j,1};
        Production{cnt,3}=str2num(MetaHIT_BileAcid_Production.(mets{i,1}){j,2});
        group=MetaHIT_IndividualIDs(find(strcmp(MetaHIT_BileAcid_Production.(mets{i,1}){j,1},MetaHIT_IndividualIDs(:,1))),2);
        if strcmp(group{1},'healthy')
            Production{cnt,2}='Healthy_MetaHIT';
        end
        if strcmp(group{1},'CD')
            Production{cnt,2}='CD_MetaHIT';
        end
        if strcmp(group{1},'UC')
            Production{cnt,2}='UC_MetaHIT';
        end
        cnt=cnt+1;
    end
    
    load('HMP_BileAcid_Production.mat');
    for j=2:size(HMP_BileAcid_Production.(mets{i,1}),1)
        Production{cnt,1}=HMP_BileAcid_Production.(mets{i,1}){j,1};
        Production{cnt,3}=str2num(HMP_BileAcid_Production.(mets{i,1}){j,2});
        Production{cnt,2}='Healthy_HMP';
        cnt=cnt+1;
    end
    % reaction abundance
    RxnAbundance={};
    RxnAbundance{1,1}='Individual';
    RxnAbundance{1,2}='Stratification';
    cnt=2;
    load(strcat(savePath,'pIBD_ReactionAbundance.mat'));
    load('pIBD_IndividualIDs.mat');
    for k=2:size(ReactionAbundance.('Total'),1)
            RxnAbundance{cnt,1}=ReactionAbundance.('Total'){k,1};
            group=pIBD_IndividualIDs(find(strcmp(ReactionAbundance.('Total'){k,1},pIBD_IndividualIDs(:,1))),2);
            if strcmp(group{1},'healthy')
                RxnAbundance{cnt,2}='Healthy_pIBD';
            end
            if strcmp(group{1},'CD')
                RxnAbundance{cnt,2}='CD_pIBD';
            end
            for j=2:size(ReactionAbundance.('Total'),2)
                RxnAbundance{1,j+1}=ReactionAbundance.('Total'){1,j};
                RxnAbundance{cnt,j+1}=ReactionAbundance.('Total'){k,j};
            end
            cnt=cnt+1;
    end
    load(strcat(savePath,'MetaHIT_ReactionAbundance.mat'));
    load('MetaHIT_BileAcid_Production.mat');
    load('MetaHIT_IndividualIDs.mat');
    % order doesn't always match
    for j=2:size(MetaHIT_BileAcid_Production.(mets{i,1}),1)
        indID=find(strcmp(MetaHIT_BileAcid_Production.(mets{i,1}){j,1},ReactionAbundance.('Total')(:,1)));
        RxnAbundance{cnt,1}=ReactionAbundance.('Total'){indID,1};
        group=MetaHIT_IndividualIDs(find(strcmp(ReactionAbundance.('Total'){indID,1},MetaHIT_IndividualIDs(:,1))),2);
        if strcmp(group{1},'healthy')
            RxnAbundance{cnt,2}='Healthy_MetaHIT';
        end
        if strcmp(group{1},'CD')
            RxnAbundance{cnt,2}='CD_MetaHIT';
        end
        if strcmp(group{1},'UC')
            RxnAbundance{cnt,2}='UC_MetaHIT';
        end
        for k=2:size(ReactionAbundance.('Total'),2)
            RxnAbundance{cnt,k+1}=ReactionAbundance.('Total'){indID,k};
        end
        cnt=cnt+1;
    end
    load(strcat(savePath,'HMP_ReactionAbundance.mat'));
    for k=2:size(ReactionAbundance.('Total'),1)
        RxnAbundance{cnt,1}=ReactionAbundance.('Total'){k,1};
        RxnAbundance{cnt,2}='Healthy_HMP';
        for j=2:size(ReactionAbundance.('Total'),2)
            RxnAbundance{cnt,j+1}=ReactionAbundance.('Total'){k,j};
        end
        cnt=cnt+1;
    end
end

load('All_ShadowPricesExtracted.mat');
% extract and summarize in table-data shown in Table S12
All_ShadowPrices_Table={};
for i=2:size(ShadowPricesExtracted,1)
    All_ShadowPrices_Table{1,i}=ShadowPricesExtracted{i,1};
end
rowCnt=2;
for j=2:size(ShadowPricesExtracted,2)
    SPmets={};
    cnt=1;
    for i=2:size(ShadowPricesExtracted,1)
        if ~isempty(ShadowPricesExtracted{i,j})
        for k=1:size(ShadowPricesExtracted{i,j},1)
            SPmets{cnt,1}=ShadowPricesExtracted{i,j}{k,1};
            cnt=cnt+1;
        end
        end
    end
    SPmets=unique(SPmets);
    for k=1:length(SPmets)
        All_ShadowPrices_Table{rowCnt,1}=strcat(SPmets{k},'_',ShadowPricesExtracted{1,j});
        for i=2:size(ShadowPricesExtracted,1)
            if ~isempty(ShadowPricesExtracted{i,j})
                findMet=find(strcmp(ShadowPricesExtracted{i,j}(:,1),SPmets{k}));
                if ~isempty(findMet)
                    All_ShadowPrices_Table{rowCnt,i}=ShadowPricesExtracted{i,j}{findMet,2};
                else
                    All_ShadowPrices_Table{rowCnt,i}=0;
                end
            else
                All_ShadowPrices_Table{rowCnt,i}=0;
            end
        end
        rowCnt=rowCnt+1;
    end
end
save(strcat(savePath,'All_ShadowPrices_Table'),'All_ShadowPrices_Table');
