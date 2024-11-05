
% create data for Figure 2

clear all
rootDir = pwd;

%% export some statistics on taxa in the resource
% get statistics for taxon presence
taxa={'Phylum','Class','Order','Family','Genus'};
taxonStats={'','Pasolli','Almeida','APOLLO'};
taxonStats{2,1}='Strains';
for i=1:length(taxa)
    taxonStats{i+2,1}=taxa{i};
end
taxonStats{8,1}='Named species';

% define the two resources
resources={
    'Pasolli'
    'Almeida'
    };

versions = {
    '_draft','translatedDraftReconstructions'
    '_refined','refinedReconstructions'
    };

mkdir([rootDir filesep 'results'])
mkdir([rootDir filesep 'results' filesep 'strains'])
mkdir([rootDir filesep 'results' filesep 'strains' filesep 'StatsByPhylum'])

% get taxon stats for all reconstructions taken together
taxonomyAll = readInputTableForPipeline([rootDir filesep 'input' filesep 'Combined_taxonomy_info.txt']);
taxonStats{2,4}=length(unique(taxonomyAll(2:end,1)));
for i=1:length(taxa)
    taxCol=find(strcmp(taxonomyAll(1,:),taxa{i}));
    getTax=unique(taxonomyAll(2:end,taxCol));
    getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
    taxonStats{i+2,4}=length(getTax);
end
% get number of characterized and uncharacterized species
taxCol=find(strcmp(taxonomyAll(1,:),'Species'));
getTax=unique(taxonomyAll(2:end,taxCol));
getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
chrTax=getTax;
chrTax(find(contains(chrTax,'sp.')),:)=[];
chrTax(find(contains(chrTax,'nov.')),:)=[];
chrTax(find(contains(chrTax,'cf')),:)=[];
chrTax(find(contains(chrTax,'iales')),:)=[];
chrTax(find(contains(chrTax,'aceae')),:)=[];
chrTax(find(contains(chrTax,'uncultured')),:)=[];
chrTax(find(contains(chrTax,'taxon')),:)=[];
chrTax(find(strncmp(chrTax,'bacterium',9)),:)=[];

taxonStats{8,4}=length(chrTax);
taxonStats{9,4}=length(getTax) - length(chrTax);
cell2csv([rootDir filesep 'results' filesep 'strains' filesep 'Taxon_statistics.csv'],taxonStats);

%% create a table with stats on both resources separately and together
table={'','Almeida','Pasolli','APOLLO','AGORA2'
    'Number of reconstructions','','','',''
    'Aerobic growth','','','',''
    'Anaerobic growth','','','',''
    'Aerobic ATP production','','','',''
    'Anaerobic ATP production','','','',''
    'Number of reactions','','','',''
    'Number of metabolites','','','',''
    'Number of genes','','','',''
    'Flux consistent reactions','','','',''
    'Stoich consistent reactions','','','',''
    'Unique reactions','','','',''
    'Unique metabolites','','','',''
    };

% load data
load([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'All_statistics_' resources{1,1} '_refined.mat']);
data_Pasolli = stats;
data_Pasolli(1,:) = [];
load([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'All_statistics_' resources{2,1} '_refined.mat']);
data_Almeida = stats;
data_Almeida(1,:) = [];
data_combined=vertcat(data_Pasolli,data_Almeida(1:end,:));
load([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'All_statistics_AGORA2.mat']);
data_AGORA2 = stats;
data_AGORA2(1,:) = [];
data_Pasolli = cell2mat(data_Pasolli(:,2:end));
data_Almeida = cell2mat(data_Almeida(:,2:end));
data_combined = cell2mat(data_combined(:,2:end));
data_AGORA2 = cell2mat(data_AGORA2(:,2:end));

table{2,2}=length(data_Almeida);
table{2,3}=length(data_Pasolli);
table{2,4}=length(data_combined);
table{2,5}=length(data_AGORA2);

% there may be NaNs in the data
data_Almeida(find(isnan(data_Almeida(:,6))),:)=[];
data_Almeida(find(isnan(data_Almeida(:,7))),:)=[];
data_Pasolli(find(isnan(data_Pasolli(:,6))),:)=[];
data_Pasolli(find(isnan(data_Pasolli(:,7))),:)=[];
data_combined(find(isnan(data_combined(:,6))),:)=[];
data_combined(find(isnan(data_combined(:,7))),:)=[];

for i=3:4
% Almeida averages
av = mean(data_Almeida(:,i+1));
s = std(data_Almeida(:,i+1));
table{i,2} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
% Pasolli averages
av = mean(data_Pasolli(:,i+1));
s = std(data_Pasolli(:,i+1));
table{i,3} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
% total averages
av = mean(data_combined(:,i+1));
s = std(data_combined(:,i+1));
table{i,4} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
% AGORA2 averages
av = mean(data_AGORA2(:,i+1));
s = std(data_AGORA2(:,i+1));
table{i,5} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
end

for i=5:11
    % Almeida averages
    av = mean(data_Almeida(:,i+4));
    s = std(data_Almeida(:,i+4));
    table{i,2} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
    % Pasolli averages
    av = mean(data_Pasolli(:,i+4));
    s = std(data_Pasolli(:,i+4));
    table{i,3} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
    % total averages
    av = mean(data_combined(:,i+4));
    s = std(data_combined(:,i+4));
    table{i,4} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
    % AGORA2 averages
    av = mean(data_AGORA2(:,i+4));
    s = std(data_AGORA2(:,i+4));
    table{i,5} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
end
cell2csv([rootDir filesep 'results' filesep 'strains' filesep 'Summary_ReconStats.csv'],table);

%% get stats on taxon distribution

taxTable={'','Almeida','Pasolli','APOLLO','AGORA2'
    'Phyla','','','',''
    'Classes','','','',''
    'Order','','','',''
    'Families','','','',''
    'Genera','','','',''
    'Species','','','',''
    'Unclassified on species level','','','',''
};
taxa={'Phylum','Class','Order','Family','Genus','Species'};

infoFiles = {
    [rootDir filesep 'input' filesep 'Pasolli_genomes_taxonomy_info.txt']
    [rootDir filesep 'input' filesep 'Almeida_genomes_taxonomy_info.txt']
    [rootDir filesep 'input' filesep 'Combined_taxonomy_info.txt']
    'AGORA2_infoFile.xlsx'
    };

for i=1:length(infoFiles)
    info = readInputTableForPipeline(infoFiles{i});
    for j=1:length(taxa)
        taxCol=find(strcmp(info(1,:),taxa{j}));
        getTax=unique(info(2:end,taxCol));
        getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
        getTax(find(strcmp(getTax,'')),:)=[];
        getTax(find(strcmp(getTax,'N/A')),:)=[];
        getTax(find(strcmp(getTax,'NA')),:)=[];
        taxTable{j+1,i+1}=length(getTax);
    end
    % get unclassified at least on species level
    taxCol=find(strcmp(info(1,:),taxa{j}));
    uncl=0;
    getUncl=find(strncmp(info(2:end,taxCol),'unclassified',length('unclassified')));
    uncl = uncl + length(getUncl);
    getUncl=find(strcmp(info(2:end,taxCol),''));
    uncl = uncl + length(getUncl);
    getUncl=find(strcmp(info(2:end,taxCol),'N/A'));
    uncl = uncl + length(getUncl);
    getUncl=find(strcmp(info(2:end,taxCol),'NA'));
    uncl = uncl + length(getUncl);
    taxTable{8,i+1} = uncl/(size(info,1)-1);
end
cell2csv([rootDir filesep 'results' filesep 'strains' filesep 'Summary_TaxStats.csv'],table);

%% separate the results by phylum for Figure 2g-h
taxonomyAll = readInputTableForPipeline([rootDir filesep 'input' filesep 'Combined_taxonomy_info.txt']);

[phyla,~,J]=unique(taxonomyAll(2:end,find(strcmp(taxonomyAll(1,:),'Phylum'))));
% remove the phyla with only few entries and unclassified bacteria
cnt = histc(J, 1:numel(phyla));
phyla(cnt<3)=[];
phyla(strncmp(phyla,'unclassified',length('unclassified')))=[];
phyla(strncmp(phyla,'Unclassified',length('Unclassified')))=[];
phyla(strncmp(phyla,'N/A',3))=[];
phyla(strcmp(phyla,''))=[];

% get the statistics by phylum
% load data
load([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'All_statistics_' resources{1,1} '_refined.mat']);
data_Pasolli = stats;
load([rootDir filesep 'data' filesep 'plot_ModelStatistics' filesep 'All_statistics_' resources{2,1} '_refined.mat']);
data_Almeida = stats;
data_combined=vertcat(data_Pasolli,data_Almeida(1:end,:));

for i=2:size(data_combined,2)
    if strcmp(data_combined{1,i},'Growth_aerobic_CM') || strcmp(data_combined{1,i},'Reactions')
        statsByPhylum={};
        for j=1:length(phyla)
            statsByPhylum{j,1}=phyla{j};
            allbacs=taxonomyAll(find(strcmp(taxonomyAll(:,find(strcmp(taxonomyAll(1,:),'Phylum'))),phyla{j})),1);
            for k=1:length(allbacs)
                statsByPhylum{j,k+1}=data_combined{find(strcmp(data_combined(:,1),allbacs{k})),i};
            end
        end
        statsByPhylum=flip(statsByPhylum);
        statsByPhylum=statsByPhylum';
        writetable(cell2table(statsByPhylum),[rootDir filesep 'results' filesep 'strains' filesep 'StatsByPhylum' filesep data_combined{1,i} '_combined.csv'],'writeVariableNames',false);
    end
end

% The plots can subsequently be generated by uploading the generated files
% at http://shiny.chemgrid.org/boxplotr/.
