
% get overlap in number of taxa between APOLLO/AGORA2

clear all
rootDir = pwd;

infoFilePasolli = readInputTableForPipeline([rootDir filesep 'input' filesep 'Pasolli_genomes_taxonomy_info.txt']);
infoFileAlmeida = readInputTableForPipeline([rootDir filesep 'input' filesep 'Almeida_genomes_taxonomy_info.txt']);
infoFileAGORA2 = readInputTableForPipeline('AGORA2_infoFile.xlsx');

taxLevels = {'Phylum','Class','Order','Family','Genus','Species'};

for i=1:length(taxLevels)
    overlap={'Almeida','Pasolli','AGORA2'};
    % Almeida
    taxCol = find(strcmp(infoFileAlmeida(1,:),taxLevels{i}));
    taxaAlmeida = unique(infoFileAlmeida(2:end,taxCol));
    taxaAlmeida(find(contains(taxaAlmeida, 'unclassified'))) = [];
    taxaAlmeida(strcmp(taxaAlmeida, '')) = [];
    overlap(2:length(taxaAlmeida)+1,1) = taxaAlmeida;
    %Pasolli
    taxCol = find(strcmp(infoFilePasolli(1,:),taxLevels{i}));
    taxaPasolli = unique(infoFilePasolli(2:end,taxCol));
    taxaPasolli(find(contains(taxaPasolli, 'unclassified'))) = [];
    taxaPasolli(strcmp(taxaPasolli, '')) = [];
    overlap(2:length(taxaPasolli)+1,2) = taxaPasolli;
    % AGORA2
    taxCol = find(strcmp(infoFileAGORA2(1,:),taxLevels{i}));
    taxaAGORA2 = unique(infoFileAGORA2(2:end,taxCol));
    taxaAGORA2(find(contains(taxaAGORA2, 'unclassified'))) = [];
    taxaAGORA2(strcmp(taxaAGORA2, '')) = [];
    overlap(2:length(taxaAGORA2)+1,3) = taxaAGORA2;
    % export overlap
    writetable(cell2table(overlap),[rootDir filesep 'results' filesep 'strains' filesep taxLevels{i} '_overlap.csv'],'writeVariableNames',false)
end

% afterwards, the plot is generated in https://bioinformatics.psb.ugent.be/webtools/Venn/

% then get number of unique taxa

infoFileAPOLLO = readInputTableForPipeline(['input' filesep 'Combined_taxonomy_info.txt']);

taxStats={'','Almeida','Pasolli','APOLLO','AGORA2','APOLLO + AGORA2'};

for i=1:length(taxLevels)
    taxStats{i+1,1} = taxLevels{i};
    % Almeida
    taxCol = find(strcmp(infoFileAlmeida(1,:),taxLevels{i}));
    taxaAlmeida = unique(infoFileAlmeida(2:end,taxCol));
    if i==6
        taxStats{i+2,1} = 'Unclassified on species level';
        taxStats{i+2,2} = length(find(strcmp(infoFileAlmeida(2:end,taxCol),'')))/(size(infoFileAlmeida,1)-1);
    end
    taxaAlmeida(find(contains(taxaAlmeida, 'unclassified'))) = [];
    taxaAlmeida(strcmp(taxaAlmeida, '')) = [];
    taxStats{i+1,2} = length(taxaAlmeida);
    %Pasolli
    taxCol = find(strcmp(infoFilePasolli(1,:),taxLevels{i}));
    taxaPasolli = unique(infoFilePasolli(2:end,taxCol));
    if i==6
        taxStats{i+2,3} = length(find(strcmp(infoFilePasolli(2:end,taxCol),'')))/(size(infoFilePasolli,1)-1);
    end
    taxaPasolli(find(contains(taxaPasolli, 'unclassified'))) = [];
    taxaPasolli(strcmp(taxaPasolli, '')) = [];
    taxStats{i+1,3} = length(taxaPasolli);
    % APOLLO
    taxCol = find(strcmp(infoFileAPOLLO(1,:),taxLevels{i}));
    taxaAPOLLO = unique(infoFileAPOLLO(2:end,taxCol));
    if i==6
        taxStats{i+2,4} = length(find(strcmp(infoFileAPOLLO(2:end,taxCol),'')))/(size(infoFileAPOLLO,1)-1);
        unTaxAPOLLO=length(find(strcmp(infoFileAPOLLO(2:end,taxCol),'')));
    end
    taxaAPOLLO(find(contains(taxaAPOLLO, 'unclassified'))) = [];
    taxaAPOLLO(strcmp(taxaAPOLLO, '')) = [];
    taxStats{i+1,4} = length(taxaAPOLLO);
    % AGORA2
    taxCol = find(strcmp(infoFileAGORA2(1,:),taxLevels{i}));
    taxaAGORA2 = unique(infoFileAGORA2(2:end,taxCol));
    if i==6
        taxStats{i+2,5} = length(find(strcmp(infoFileAGORA2(2:end,taxCol),'')))/(size(infoFileAGORA2,1)-1);
        unTaxAGORA2=length(find(strcmp(infoFileAGORA2(2:end,taxCol),'')));
    end
    taxaAGORA2(find(contains(taxaAGORA2, 'unclassified'))) = [];
    taxaAGORA2(strcmp(taxaAGORA2, '')) = [];
    taxStats{i+1,5} = length(taxaAGORA2);
    % APOLLO + AGORA2
    taxaAll = union(taxaAGORA2, taxaAPOLLO);
    taxStats{i+1,6} = length(taxaAll);
    if i==6
        taxStats{i+2,6} = (unTaxAPOLLO+unTaxAGORA2)/((size(infoFileAPOLLO,1)-1)+(size(infoFileAGORA2,1)-1));
    end
end
