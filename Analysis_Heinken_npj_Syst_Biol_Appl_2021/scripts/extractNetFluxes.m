% This script extracts the computed net uptake and secretion fluxes for the
% personalized microbiome models generated for the PLEASE/COMBO cohort.
% Almut Heinken, 03/21

metadata=table2cell(readtable(metadataPath, 'ReadVariableNames', false));

% first get the list of all possible exchange reactions and make a table
DietExchanges={};
for i=2:size(metadata,1)
    i
    % load model
     load([modelPath filesep metadata{i,1} '.mat']);
    DietExRxns=model.rxns(find(strncmp('Diet_EX_',model.rxns,8)));
    DietExchanges=vertcat(DietExchanges,DietExRxns);
end
DietExchanges=unique(DietExchanges);
%% Calculate statistics on uptake/secretion
% First extract absolute fluxes into tables
for i=2:size(metadata,1)
    i
    % load model
     load([modelPath filesep metadata{i,1} '.mat']);
   % load corresponding uptake and secretion fluxes
    load([fluxPath filesep 'summary_' metadata{i,1} '.mat']);
    DietExchangeMinFluxes{i,1}=metadata{i,1};
    DietExchangeMaxFluxes{i,1}=metadata{i,1};
    FecalExchangeMinFluxes{i,1}=metadata{i,1};
    FecalExchangeMaxFluxes{i,1}=metadata{i,1};
    for k=1:length(DietExchanges)
        % biomass causes problems
        if ~strcmp(DietExchanges{k},'EX_biomass[c]') && ~strcmp(DietExchanges{k},'Diet_EX_biomass[c]')
            currExchDiet=find(strcmp(model.rxns,DietExchanges{k}));
            FecalExch=regexprep(DietExchanges{k},'Diet_EX_','EX_');
            FecalExch=regexprep(FecalExch,'\[d]','\[fe]');
            currExchFecal=find(strcmp(model.rxns,FecalExch));
            DietExchangeMinFluxes{1,k+1}=regexprep(DietExchanges{k},'Diet_EX_','');
            DietExchangeMinFluxes{1,k+1}=regexprep(DietExchangeMinFluxes{1,k+1},'\[d]','');
            DietExchangeMaxFluxes{1,k+1}=DietExchangeMinFluxes{1,k+1};
            FecalExchangeMinFluxes{1,k+1}=DietExchangeMinFluxes{1,k+1};
            FecalExchangeMaxFluxes{1,k+1}=DietExchangeMinFluxes{1,k+1};
            if ~isempty(currExchDiet)
                DietExchangeMinFluxes{i,k+1}=(sumResults(find(sumResults(:,1)==currExchDiet),2));
                DietExchangeMaxFluxes{i,k+1}=(sumResults(find(sumResults(:,1)==currExchDiet),3));
                FecalExchangeMinFluxes{i,k+1}=(sumResults(find(sumResults(:,1)==currExchFecal),2));
                FecalExchangeMaxFluxes{i,k+1}=(sumResults(find(sumResults(:,1)==currExchFecal),3));
            else
                DietExchangeMinFluxes{i,k+1}=0;
                DietExchangeMaxFluxes{i,k+1}=0;
                FecalExchangeMinFluxes{i,k+1}=0;
                FecalExchangeMaxFluxes{i,k+1}=0;
            end
        end
    end
end
% remove empty columns
findEmptyCol=find(cellfun(@isempty,DietExchangeMinFluxes(1,:)));
DietExchangeMinFluxes(:,findEmptyCol(2:end))=[];
DietExchangeMaxFluxes(:,findEmptyCol(2:end))=[];
FecalExchangeMinFluxes(:,findEmptyCol(2:end))=[];
FecalExchangeMaxFluxes(:,findEmptyCol(2:end))=[];

save([netFluxesPath filesep 'DietExchangeMinFluxes.mat'],'DietExchangeMinFluxes');
save([netFluxesPath filesep 'DietExchangeMaxFluxes.mat'],'DietExchangeMaxFluxes');
save([netFluxesPath filesep 'FecalExchangeMinFluxes.mat'],'FecalExchangeMinFluxes');
save([netFluxesPath filesep 'FecalExchangeMaxFluxes.mat'],'FecalExchangeMaxFluxes');

% Then calculate the uptake/secretion statistics
MetaboliteStatistics{1,2}='Only_uptake';
MetaboliteStatistics{1,3}='Only_secretion';
MetaboliteStatistics{1,4}='Uptake_and_secretion';
MetaboliteStatistics{1,5}='Total_production';
MetaboliteStatistics{1,6}='No_flux';
for i=2:size(DietExchangeMinFluxes,2)
    MetaboliteStatistics(i,1)=metaboliteInfo(find(strcmp(DietExchangeMinFluxes{1,i},metaboliteInfo(:,2))),1);
    MetaboliteStatistics{i,2}=0;
    MetaboliteStatistics{i,3}=0;
    MetaboliteStatistics{i,4}=0;
    MetaboliteStatistics{i,5}=0;
    MetaboliteStatistics{i,6}=0;
    for j=2:size(DietExchangeMinFluxes,1)
        if DietExchangeMinFluxes{j,i}<-0.0000001 && abs(FecalExchangeMaxFluxes{j,i}) <= abs(DietExchangeMinFluxes{j,i})
            MetaboliteStatistics{i,2}=MetaboliteStatistics{i,2}+1;
        elseif abs(DietExchangeMinFluxes{j,i})<0.0000001 && FecalExchangeMaxFluxes{j,i}>0.0000001
            MetaboliteStatistics{i,3}=MetaboliteStatistics{i,3}+1;
        elseif DietExchangeMinFluxes{j,i}<-0.0000001 && FecalExchangeMaxFluxes{j,i}>0.0000001
            MetaboliteStatistics{i,4}=MetaboliteStatistics{i,4}+1;
        elseif abs(DietExchangeMinFluxes{j,i})<0.0000001 && abs(FecalExchangeMaxFluxes{j,i})<0.0000001
            MetaboliteStatistics{i,6}=MetaboliteStatistics{i,6}+1;
        end
        MetaboliteStatistics{i,5}=MetaboliteStatistics{i,3}+MetaboliteStatistics{i,4};
    end
end
cell2csv([netFluxesPath filesep 'MetaboliteStatistics.csv'],cell(MetaboliteStatistics));
% Summarize the statistics
StatisticsSummarized{1,2}='Number';
StatisticsSummarized{2,1}='Metabolites_only_uptake';
StatisticsSummarized{3,1}='Metabolites_only_secretion';
StatisticsSummarized{4,1}='Metabolites_uptake_and_secretion';
StatisticsSummarized{5,1}='Metabolites_no_flux';
StatisticsSummarized{6,1}='Metabolites_produced_=>90';
StatisticsSummarized{7,1}='Metabolites_produced_=>25_<90';
StatisticsSummarized{8,1}='Metabolites_produced_<25';
for i=2:size(StatisticsSummarized,1)
    StatisticsSummarized{i,2}=0;
end
for i=2:size(MetaboliteStatistics,1)
    if MetaboliteStatistics{i,2}>0 && MetaboliteStatistics{i,4}==0
        StatisticsSummarized{2,2}=StatisticsSummarized{2,2}+1;
    end
    if MetaboliteStatistics{i,3}>0 && MetaboliteStatistics{i,4}==0
        StatisticsSummarized{3,2}=StatisticsSummarized{3,2}+1;
    end
    if MetaboliteStatistics{i,4}>0
        StatisticsSummarized{4,2}=StatisticsSummarized{4,2}+1;
    end
    if  MetaboliteStatistics{i,2}==0 && MetaboliteStatistics{i,3}==0 && MetaboliteStatistics{i,4}==0
        StatisticsSummarized{5,2}=StatisticsSummarized{5,2}+1;
    end
    if MetaboliteStatistics{i,5} >= max(cell2mat(MetaboliteStatistics(2:end,5)))*0.9
        StatisticsSummarized{6,2}=StatisticsSummarized{6,2}+1;
    end
    if MetaboliteStatistics{i,5} >= max(cell2mat(MetaboliteStatistics(2:end,5)))*0.25 &&  MetaboliteStatistics{i,5} < max(cell2mat(MetaboliteStatistics(2:end,5)))*0.9
        StatisticsSummarized{7,2}=StatisticsSummarized{7,2}+1;
    end
    if MetaboliteStatistics{i,5} >0 && MetaboliteStatistics{i,5} < max(cell2mat(MetaboliteStatistics(2:end,5)))*0.25
        StatisticsSummarized{8,2}=StatisticsSummarized{8,2}+1;
    end
end
cell2csv([netFluxesPath filesep 'StatisticsSummarized.csv'],cell(StatisticsSummarized));
%% Calculate net uptake and net production potential of each microbiome
% Extract net fluxes into table

netProd={};
netUpt={};
for i=2:size(metadata,1)
    i
   % load model
   load([modelPath filesep metadata{i,1} '.mat']);
   % load corresponding uptake and secretion fluxes
    load([fluxPath filesep 'summary_' metadata{i,1} '.mat']);
    netProd{i,1}=metadata{i,1};
    for k=1:length(DietExchanges)
        % biomass causes problems
        if ~strcmp(DietExchanges{k},'EX_biomass[c]') && ~strcmp(DietExchanges{k},'Diet_EX_biomass[c]')
            currExchDiet=find(strcmp(model.rxns,DietExchanges{k}));
            FecalExch=regexprep(DietExchanges{k},'Diet_EX_','EX_');
            FecalExch=regexprep(FecalExch,'\[d]','\[fe]');
            currExchFecal=find(strcmp(model.rxns,FecalExch));
            netProd{1,k+1}=FecalExch;
            if ~isempty(currExchDiet)
                netProd{i,k+1}=(sumResults(find(sumResults(:,1)==currExchFecal),3))+(sumResults(find(sumResults(:,1)==currExchDiet),2));
                if  netProd{i,k+1}<0
                    netProd{i,k+1}=0;
                end
            else
                netProd{i,k+1}=0;
            end
        end
    end
    netUpt{i,1}=metadata{i,1};
    for k=1:length(DietExchanges)
        % biomass causes problems
        if ~strcmp(DietExchanges{k},'EX_biomass[c]') && ~strcmp(DietExchanges{k},'Diet_EX_biomass[c]')
            currExchDiet=find(strcmp(model.rxns,DietExchanges{k}));
            FecalExch=regexprep(DietExchanges{k},'Diet_EX_','EX_');
            FecalExch=regexprep(FecalExch,'\[d]','\[fe]');
            currExchFecal=find(strcmp(model.rxns,FecalExch));
            netUpt{1,k+1}=FecalExch;
            if ~isempty(currExchDiet)
                netUpt{i,k+1}=abs((sumResults(find(sumResults(:,1)==currExchFecal),2))+(sumResults(find(sumResults(:,1)==currExchDiet),3)));
            else
                netUpt{i,k+1}=0;
            end
        end
    end
end
cell2csv([netFluxesPath filesep 'netProductionFluxes.csv'],cell(netProd'));
cell2csv([netFluxesPath filesep 'netUptakeFluxes.csv'],cell(netUpt'));

%% export the fluxes
delInd=find(cellfun(@isempty,netProd(1,:)));
netProd{1,1}='Metabolite';
netProd(:,delInd(2:end))=[];
netProd(1,:) = strrep(netProd(1,:),'EX_','');
netProd(1,:) = strrep(netProd(1,:),'[fe]','');

% delete zero columns
delArray=[];
cnt=1;
for i=2:size(netProd,2)
    if abs(sum(cell2mat(netProd(2:end,i)))) <0.000001
        delArray(cnt,1)=i;
        cnt=cnt+1;
    end
end
netProd(:,delArray)=[];
for i=2:size(netProd,2)
    netProd(1,i)=metaboliteInfo(find(strcmp(metaboliteInfo(:,2),netProd{1,i})),1);
end
cell2csv([netFluxesPath filesep 'netProductionFluxes_Exported.csv'],cell(netProd'));

delInd=find(cellfun(@isempty,netUpt(1,:)));
netUpt{1,1}='Metabolite';
netUpt(:,delInd(2:end))=[];
netUpt(1,:) = strrep(netUpt(1,:),'EX_','');
netUpt(1,:) = strrep(netUpt(1,:),'[fe]','');

% delete zero columns and values that are too similar
delArray=[];
cnt=1;
for i=2:size(netUpt,2)
    if abs(sum(cell2mat(netUpt(2:end,i)))) <0.000001 || std(cell2mat(netUpt(2:end,i)))<0.000001
        delArray(cnt,1)=i;
        cnt=cnt+1;
    end
end
netUpt(:,delArray)=[];
for i=2:size(netUpt,2)
    netUpt(1,i)=metaboliteInfo(find(strcmp(metaboliteInfo(:,2),netUpt{1,i})),1);
end
cell2csv([netFluxesPath filesep 'netUptakeFluxes_Exported.csv'],cell(netUpt'));
