
% map taxa found in Table S9 of Yachida et al 2019 (PMID:31171880) to the
% corresponding AGORA2 pan-species model IDs.

% define taxa with updated nomenclature that need to be translated
speciesMapped={
    'Actinomyces_cardiffensis','Schaalia_cardiffensis'
    'Actinomyces_europaeus','Gleimia_europaea'
    'Actinomyces_odontolyticus','Schaalia_odontolytica'
    'Actinomyces_turicensis','Schaalia_turicensis'
    'Alistipes_sp_AP11','Alistipes_ihumii'
    'Arcobacter_butzleri','Aliarcobacter_butzleri'
    'Atopobium_parvulum','Lancefieldella_parvula'
    'Atopobium_rimae','Lancefieldella_rimae'
    'Bacteroidales_bacterium_ph8','Alistipes_communis'
    'Bacteroides_barnesiae','Phocaeicola_barnesiae'
    'Bacteroides_coprocola','Phocaeicola_coprocola'
    'Bacteroides_coprophilus','Phocaeicola_coprophilus'
    'Bacteroides_dorei','Phocaeicola_dorei'
    'Bacteroides_massiliensis','Phocaeicola_massiliensis'
    'Bacteroides_plebeius','Phocaeicola_plebeius'
    'Bacteroides_salanitronis','Phocaeicola_salanitronis'
    'Clostridium_asparagiforme','Enterocloster_asparagiformis'
    'Clostridium_bartlettii','Intestinibacter_bartlettii'
    'Clostridium_bifermentans','Paraclostridium_bifermentans'
    'Clostridium_bolteae','Enterocloster_bolteae'
    'Clostridium_citroniae','Enterocloster_citroniae'
    'Clostridium_clostridioforme','Enterocloster_clostridioformis'
    'Clostridium_difficile','Clostridioides_difficile'
    'Clostridium_glycolicum','Terrisporobacter_glycolicus'
    'Clostridium_hathewayi','Hungatella_hathewayi'
    'Clostridium_hiranonis','Peptacetobacter_hiranonis'
    'Clostridium_innocuum','Erysipelatoclostridium_innocuum'
    'Clostridium_nexile','Tyzzerella_nexilis'
    'Clostridium_ramosum','Erysipelatoclostridium_ramosum'
    'Clostridium_sordellii','Paeniclostridium_sordellii'
    'Coprobacillus_sp_29_1','Coprobacillus_cateniformis'
    'Coriobacteriaceae_bacterium_phI','Enorma_massiliensis'
    'Eggerthella_sp_HGA1','Eggerthella_lenta'
    'Enterobacter_aerogenes','Klebsiella_aerogenes'
    'Enterorhabdus_caecimuris','Adlercreutzia_caecimuris'
    'Escherichia_hermannii','Atlantibacter_hermannii'
    'Eubacterium_biforme','Holdemanella_biformis'
    'Eubacterium_cylindroides','Faecalitalea_cylindroides'
    'Eubacterium_dolichum','Amedibacillus_dolichus'
    'Eubacterium_eligens','Lachnospira_eligens'
    'Eubacterium_hallii','Anaerobutyricum_hallii'
    'Holdemania_sp_AP2','Holdemania_massiliensis'
    'Lachnospiraceae_oral_taxon_107','Lachnospiraceae_bacterium_oral_taxon_500'
    'Lactobacillus_brevis','Levilactobacillus_brevis'
    'Lactobacillus_curvatus','Latilactobacillus_curvatus'
    'Lactobacillus_fermentum','Limosilactobacillus_fermentum'
    'Lactobacillus_gastricus','Limosilactobacillus_gastricus'
    'Lactobacillus_mucosae','Limosilactobacillus_mucosae'
    'Lactobacillus_oris','Limosilactobacillus_oris'
    'Lactobacillus_pentosus','Lactiplantibacillus_pentosus'
    'Lactobacillus_plantarum','Lactiplantibacillus_plantarum'
    'Lactobacillus_reuteri','Limosilactobacillus_reuteri'
    'Lactobacillus_rhamnosus','Lacticaseibacillus_rhamnosus'
    'Lactobacillus_ruminis','Ligilactobacillus_ruminis'
    'Lactobacillus_sakei','Latilactobacillus_sakei'
    'Lactobacillus_salivarius','Ligilactobacillus_salivarius'
    'Lactobacillus_sanfranciscensis','Fructilactobacillus_sanfranciscensis'
    'Lactobacillus_vaginalis','Limosilactobacillus_vaginalis'
    'Leuconostoc_gasicomitatum','Leuconostoc_gelidum'
    'Prevotella_sp_oral_taxon_473','Alloprevotella_sp_oral_taxon_473'
    'Propionibacterium_acnes','Cutibacterium_acnes'
    'Propionibacterium_propionicum','Arachnia_propionica'
    'Ruminococcus_obeum','Blautia_obeum'
    'Streptococcus_oligofermentans','Streptococcus_cristatus'
    'Streptococcus_tigurinus','Streptococcus_oralis'
    };

coverage = readInputTableForPipeline([rootDir filesep 'input' filesep 'Table_S9_Yachida_2019.csv']);
coverage(:,1) = strrep(coverage(:,1),' ','_');

for i=1:length(speciesMapped)
    coverage(:,1) = strrep(coverage(:,1),speciesMapped{i,1},speciesMapped{i,2});
end

infoFile = readInputTableForPipeline('AGORA2_infoFile.xlsx');
species=unique(infoFile(2:end,5));
species = strrep(species,' ','_');
species = strrep(species,'.','');
species = strrep(species,'-','_');
species = strrep(species,'/','_');
species = strrep(species,'(','');
species = strrep(species,')','');
species = strrep(species,'[','');
species = strrep(species,']','');

% remove taxa not in AGORA2
[C,IA]=setdiff(coverage(:,1),species,'stable');
coverage(IA(2:end),:)=[];

% adjust namespaces
for i=2:size(coverage,1)
    coverage{i,1} = ['pan' coverage{i,1}];
end
for i=2:size(coverage,2)
    coverage{1,i} = ['sample' coverage{1,i}];
end
coverage(1,2:end)=strrep(coverage(1,2:end),'.','_');
cell2csv([rootDir filesep 'input' filesep 'mappedCoverage_CRC.csv'],coverage)

% normalize coverage
abunFilePath = [rootDir filesep 'input' filesep 'mappedCoverage_CRC.csv'];
[normalizedCoverage,normalizedCoveragePath] = normalizeCoverage(abunFilePath,0);
