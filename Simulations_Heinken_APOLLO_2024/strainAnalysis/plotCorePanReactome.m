

clear all
loadedData = parquetread(['Model_properties_analysis' filesep '150k_90k_parquet_files' filesep 'reactionPresence_combined_refined.parquet']);

% load all strain taxonomy data
info_Pasolli=readInputTableForPipeline(['input' filesep 'Pasolli_genomes_taxonomy_info.xlsx']);
strain_names = cellstr(table2cell(loadedData(:,end)));
loadedData(:,end)=[];
rxnNames = loadedData.Properties.VariableNames;

% remove any AGORA2 strains
[~,IA]=setdiff(info_Pasolli(:,1),strain_names,'stable');
info_Pasolli(IA(2:end),:)=[];
info_Almeida=readInputTableForPipeline(['input' filesep 'Almeida_genomes_taxonomy_info.xlsx']);
% remove any AGORA2 strains
[~,IA]=setdiff(info_Almeida(:,1),strain_names,'stable');
info_Almeida(IA(2:end),:)=[];

% taxonLevels={'Species','Genus','Family','Order','Class','Phylum'};
taxonLevels={'Species'};

Core_pan_reactions=struct;

for i=1:length(taxonLevels)
    % create a table with information on pan and core reactome for each taxon
    taxStats={'Taxon','Number of taxa','Reactions in core reactome','Reactions in pan-reactome'};
    % get the list of all taxa on this taxon level
    taxCol_150k=find(strcmp(info_Pasolli(1,:),taxonLevels{i}));
    allTax=unique(info_Pasolli(2:end,taxCol_150k));
    taxCol_90k=find(strcmp(info_Almeida(1,:),taxonLevels{i}));
    allTax=union(allTax,unique(info_Almeida(2:end,taxCol_90k)));
    
    % removed unspecified taxa
    TF = endsWith(allTax,' sp');
    allTax(find(TF==true),:)=[];
    TF = endsWith(allTax,' unclassified');
    allTax(find(TF==true),:)=[];
        TF = endsWith(allTax,' bacterium');
    allTax(find(TF==true),:)=[];
    allTax(find(strcmp(allTax,'')),:)=[];
    allTax(find(strcmp(allTax,'N/A')),:)=[];
    allTax(find(strcmp(allTax,'NA')),:)=[];
    
    for j=1:length(allTax)
        j
        taxStats{j+1,1}=allTax{j};
        % get the reaction patterns for each strain in the taxon
        reacPat={};
        findIn150k=find(strcmp(info_Pasolli(:,taxCol_150k),allTax{j}));
        findIn90k=find(strcmp(info_Almeida(:,taxCol_90k),allTax{j}));
        
        % report total number of strains
        taxStats{j+1,2}=length(findIn150k)+length(findIn90k);
        
        for k=1:length(findIn150k)
            % get reactions in this strain
            findInData=find(strcmp(strain_names,info_Pasolli{findIn150k(k),1}));
            reacPat{length(reacPat)+1,1}=rxnNames(find(cell2mat(table2cell(loadedData(findInData,:)))==1));
        end
        for k=1:length(findIn90k)
            % get reactions in this strain
            findInData=find(strcmp(strain_names,info_Almeida{findIn90k(k),1}));
            reacPat{length(reacPat)+1,1}=rxnNames(find(cell2mat(table2cell(loadedData(findInData,:)))==1));
        end
        % get the pan-reactome of all strains
        pan_rxns={};
        for k=1:length(reacPat)
            pan_rxns=union(pan_rxns,reacPat{k,1});
        end
        taxStats{j+1,4}=length(pan_rxns);
        
        % get the core reactome of all strains
        core_rxns={};
        for k=1:length(pan_rxns)
            incore=1;
            for l=length(reacPat)
                if isempty(intersect(pan_rxns{k},reacPat{l}))
                    incore=0;
                end
            end
            if incore==1
                core_rxns{length(core_rxns)+1,1}=pan_rxns{k};
            end
        end
        taxStats{j+1,3}=length(core_rxns);
        fieldname=strrep(allTax{j},' ','_');
        fieldname=strrep(fieldname,'.','');
        fieldname=strrep(fieldname,'-','_');
        Core_pan_reactions.(taxonLevels{i}).(fieldname).('pan_rxns')=pan_rxns;
        Core_pan_reactions.(taxonLevels{i}).(fieldname).('core_rxns')=core_rxns;
    end

    % export the table
    writetable(cell2table(taxStats),['Model_properties_analysis' filesep 'Core_pan_reactome_' taxonLevels{i} '.csv'],'WriteVariableNames',false)
end

% export the lists of reactions
save(['Model_properties_analysis' filesep 'Core_pan_reactions.mat'],'Core_pan_reactions')

% plot the core and pan-reactome for the most abundant taxa
for i=1:length(taxonLevels)
    data=readInputTableForPipeline(['Model_properties_analysis' filesep 'Core_pan_reactome_' taxonLevels{i} '.csv']);
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
    title(['Pan and core reactome on ' lower(taxonLevels{i}) ' level'])  
    print(['Model_properties_analysis' filesep 'Pan_core_rxns_' taxonLevels{i}],'-dpng','-r300')
end
