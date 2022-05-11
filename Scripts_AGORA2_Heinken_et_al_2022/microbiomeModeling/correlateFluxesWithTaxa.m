
abundancePath = [rootDir filesep 'Modeling_CRC' filesep 'normCoverage_CRC.csv'];
fluxPath = [rootDir filesep 'Modeling_CRC' filesep 'Solutions_ShadowPrices' filesep 'AGORA2_CRC_Objectives.txt'];

taxonomy = readtable('AGORA2_infoFile.xlsx', 'ReadVariableNames', false);
taxonomy = table2cell(taxonomy);

%% calculate correlations between fluxes and taxon abundance-Japanese diet
[FluxCorrelations, PValues] = correlateFluxWithTaxonAbundanceCRC(abundancePath, fluxPath, 'Spearman');
save([rootDir filesep 'Modeling_CRC' filesep 'Correlations' filesep 'Taxon_Correlations_JD'],'FluxCorrelations');
close all

%% Get a printout for R plot

fluxLabels={'EX_sn38[fe]','Irinotecan deglucuronidated (sn38)';'EX_r406[fe]','R406 (Fostamatinib product)';'EX_5fura[fe]','5-Fluorouracil';'EX_dh5fura[fe]','5,6-Dihydro-5-Fluorouracil (5-fluorouracil product)';'EX_dfduri[fe]','2,2-Difluorodeoxyuridine (Gemcitabine product)';'EX_dihydro_digoxin[fe]','Dihydrodigoxin (Digoxin product)';'EX_ac5asa[fe]','N-Acetyl-5-aminosalicylic acid (Mesalamine product)';'EX_5asa[fe]','Mesalamine (active product of balsalazide)';'EX_nchlphncl[fe]','Nitrosochloramphenicol (Chloramphenicol product)';'EX_bvu[fe]',' (E)-5-(2-Bromovinyl)Uracil (Sorivudine product)';'EX_dopa[fe]','Dopamine';'EX_mtym[fe]','Tyramine';'EX_pcresol[fe]','p-cresol'};
toPlot={'';'EX_sn38[fe]';'EX_r406[fe]';'EX_5fura[fe]';'EX_dh5fura[fe]';'EX_dfduri[fe]';'EX_dihydro_digoxin[fe]';'EX_ac5asa[fe]';'EX_5asa[fe]';'EX_nchlphncl[fe]';'EX_bvu[fe]';'EX_dopa[fe]';'EX_tym[fe]';'EX_pcresol[fe]'};

levels=fieldnames(FluxCorrelations);
for i=1:length(levels)
    if ~strcmp(levels{i},'Phylum')
        correlations=FluxCorrelations.(levels{i});
        correlations{1,1}='';
        [~,IA] = setdiff(correlations(1,:),toPlot);
        correlations(:,IA)=[];
        for j=2:size(correlations,2)
            correlations{1,j}=fluxLabels{find(strcmp(fluxLabels,correlations{1,j})),2};
        end
        correlations(1,:)=strrep(correlations(1,:),',',' ');
        % remove correlations that are small
        cnt=1;
        del=[];
        for j=2:size(correlations,1)
            if ~any(abs(cell2mat(correlations(j,2:end)))>0.3)
                del(cnt,1)=j;
                cnt=cnt+1;
            end
        end
        correlations(del,:)=[];
        cell2csv([rootDir filesep 'Modeling_CRC' filesep 'Correlations' filesep 'JD_',levels{i},'.csv'],correlations);
        % create a taxonomy file for R plot
        createTaxonomy=taxonomy;
        % remove columns not needed
        findTaxonLevel=find(strcmp(createTaxonomy(1,:),levels{i}));
        createTaxonomy(:,1:findTaxonLevel-1)=[];
        % reduce to only unique entries
        [C,IA,IC] = unique(createTaxonomy(:,1),'stable');
        createTaxonomy=createTaxonomy(IA,1:size(createTaxonomy,2));
        % find entries not in data
        [C,IA] = setdiff(createTaxonomy(:,1),correlations(2:end,1),'stable');
        createTaxonomy(IA(2:end),:)=[];
        createTaxonomy(:,10:end)=[];
        cell2csv([rootDir filesep 'Modeling_CRC' filesep 'Correlations' filesep 'Taxonomy_JD_',levels{i},'.csv'],createTaxonomy);
    end
end

