% Get statistics on the number in taxa in AGORA 1.03 and AGORA2

stats={'','AGORA 1.03','AGORA2'
    'Strains','',''};
taxa={'Phylum','Class','Order','Family','Genus'};

resultsFolder = [rootDir filesep 'modelProperties'];

% AGORA 1.03
taxonomy = readInputTableForPipeline('AGORA_infoFile.xlsx');
stats{2,2}=length(unique(taxonomy(2:end,1)));
for i=1:length(taxa)
    taxCol=find(strcmp(taxonomy(1,:),taxa{i}));
    stats{i+2,1}=taxa{i};
    getTax=unique(taxonomy(2:end,taxCol));
    getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
    stats{i+2,2}=length(getTax);
end
% get number of characterized and uncharacterized species
taxCol=find(strcmp(taxonomy(1,:),'Species'));
getTax=unique(taxonomy(2:end,taxCol));
getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
chrTax=getTax;
chrTax(find(contains(chrTax,'sp.')),:)=[];
chrTax(find(contains(chrTax,'iales')),:)=[];
chrTax(find(contains(chrTax,'aceae')),:)=[];
stats{8,1}='Species, characterized';
stats{8,2}=length(chrTax);
stats{9,1}='Species, uncharacterized';
stats{9,2}=length(getTax) - length(chrTax);

% AGORA2
taxonomy = readInputTableForPipeline('AGORA2_infoFile.xlsx');
stats{2,3}=length(unique(taxonomy(2:end,1)));
for i=1:length(taxa)
    taxCol=find(strcmp(taxonomy(1,:),taxa{i}));
    getTax=unique(taxonomy(2:end,taxCol));
    getTax(find(strncmp(getTax,'unclassified',length('unclassified'))),:)=[];
    stats{i+2,3}=length(getTax);
end
% get number of characterized and uncharacterized species
taxCol=find(strcmp(taxonomy(1,:),'Species'));
getTax=unique(taxonomy(2:end,taxCol));
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

stats{8,3}=length(chrTax);
stats{9,3}=length(getTax) - length(chrTax);

% get number of human and mouse strains
taxCol=find(strcmp(taxonomy(1,:),'Host'));
humanStrains=length(find(contains(taxonomy(:,taxCol),'Human')));
mouseStrains=length(find(contains(taxonomy(:,taxCol),'Mouse')));

% get numbers for reference categories
rfCol=find(strcmp(taxonomy(1,:),'Reference category'));
[rc, ~, J]=unique(taxonomy(2:end,rfCol));
cnt = histc(J, 1:numel(rc));

save([resultsFolder filesep 'Taxon_statistics'],'stats')
