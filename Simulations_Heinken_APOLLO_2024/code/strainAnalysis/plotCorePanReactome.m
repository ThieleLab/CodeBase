
% plot the core and pan-reactome for the most abundant taxa for Figure S4
data=readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'Core_pan_reactome_Species.csv']);
data(1,:)=[];
% get 30 most abundant or otherwise all taxa
if size(data,1) >30
    cutoff=30;
else
    cutoff=size(data,1);
end
[A,I]=sort(cell2mat(data(:,2)),'descend');
core=cell2mat(data(I(1:cutoff),3));
pan=cell2mat(data(I(1:cutoff),4));
labels=data(I(1:cutoff),1);

for j=1:cutoff
    labels{j}=[data{I(j),1} ' (' num2str(data{I(j),2}) ')'];
end
f=figure;
bar(pan)
hold on
bar(core)
set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels);
set(gca,'TickLabelInterpreter','none');
xtickangle(45)
legend({'Pan reactome','Core reactome'},'Location','Northeast')
title(['Pan and core reactome on ' lower('Species') ' level'])
print([rootDir filesep 'results' filesep 'strains' filesep 'Pan_core_rxns_Species'],'-dpng','-r300')
