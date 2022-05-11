% test the coverage of AGORA1 and AGORA2 of colorectal cancer microbiomes
% and others (HMP, IBD)


% get overlap with species in AGORA1
taxonomy = readtable('AGORA_infoFile.xlsx', 'ReadVariableNames', false);
taxonomy = table2cell(taxonomy);
taxCol=find(strcmp(taxonomy(1,:),'Species'));
species=unique(taxonomy(2:end,taxCol));
species=strrep(species, ' ', '_');
species=strrep(species,'[','');
species=strrep(species,']','');
species=strrep(species,'/','_');
species=strrep(species,'-','_');
agoraSpecies1=strrep(species,'.','');

% get overlap with species in AGORA2
taxonomy = readtable('AGORA2_infoFile.xlsx', 'ReadVariableNames', false);
taxonomy = table2cell(taxonomy);
taxCol=find(strcmp(taxonomy(1,:),'Species'));
species=unique(taxonomy(2:end,taxCol));
species=strrep(species, ' ', '_');
species=strrep(species,'[','');
species=strrep(species,']','');
species=strrep(species,'/','_');
species=strrep(species,'-','_');
agoraSpecies2=strrep(species,'.','');

calculatedCoverage={'Cohort','AGORA 1.3 coverage (% taxa)','AGORA 2.0 coverage (% taxa)','AGORA 1.3 coverage (% abundance)','AGORA 2.0 coverage (% abundance)'};

%% get coverage for CRC models
abundance =readtable([rootDir filesep 'Modeling_CRC' filesep 'coverage_CRC.csv']);
abundance=table2cell(abundance);
% Final number of species not captured
% Taxa that are not considered: unclassified, viruses, eukaryotes
notConsidered={'Acidaminococcus_unclassified';'Acinetobacter_unclassified';'Actinobacillus_unclassified';'Aeromonas_phage_phiO18P';'Aeromonas_unclassified';'Aggregatibacter_unclassified';'Alistipes_unclassified';'Alloprevotella_unclassified';'Anaerostipes_unclassified';'Anaerotruncus_unclassified';'Arcobacter_unclassified';'Bacteroides_phage_B124_14';'Beet_mild_curly_top_virus';'Bilophila_unclassified';'Bombyx_mori_nucleopolyhedrovirus';'Brachyspira_unclassified';'Brevibacterium_unclassified';'Burkholderia_unclassified';'Butyrivibrio_unclassified';'C2likevirus_unclassified';'Candida_glabrata';'Capnocytophaga_unclassified';'Cellulophaga_unclassified';'Chicken_anemia_virus';'Citrobacter_unclassified';'Clostridium_phage_PhiS63';'Collinsella_unclassified';'Comamonas_unclassified';'Coprobacillus_unclassified';'Deinococcus_unclassified';'Dorea_unclassified';'Dysgonomonas_unclassified';'Eggerthella_unclassified';'Enterobacteria_phage_HK544';'Enterobacteria_phage_HK633';'Enterobacteria_phage_PsP3';'Enterobacteria_phage_ST104';'Enterobacteria_phage_mEp237';'Epsilon15likevirus_unclassified';'Escherichia_unclassified';'Facklamia_unclassified';'Gemella_unclassified';'Giardia_intestinalis';'Granulicatella_unclassified';'Granulicella_unclassified';'Holdemania_unclassified';'Human_herpesvirus_4';'JC_polyomavirus';'Klebsiella_phage_phiKO2';'Klebsiella_unclassified';'Lactobacillus_casei_paracasei';'Lactobacillus_phage_Lc_Nu';'Lactobacillus_phage_PL_1';'Leptotrichia_unclassified';'Leuconostoc_unclassified';'Lymphocryptovirus_unclassified';'Megamonas_unclassified';'Megasphaera_unclassified';'Methanobrevibacter_unclassified';'Mitsuokella_unclassified';'Mobiluncus_unclassified';'Mulikevirus_unclassified';'Naumovozyma_unclassified';'Neisseria_unclassified';'Odoribacter_unclassified';'Olsenella_unclassified';'Oscillibacter_unclassified';'P22likevirus_unclassified';'Pantoea_unclassified';'Parabacteroides_unclassified';'Paraprevotella_unclassified';'Parvimonas_unclassified';'Pediococcus_unclassified';'Peptostreptococcaceae_noname_unclassified';'Peptostreptococcus_unclassified';'Porcine_type_C_oncovirus';'Propionibacterium_phage_P100D';'Propionibacterium_phage_P100_A';'Propionibacterium_phage_P104A';'Propionibacterium_phage_PAD20';'Propionibacterium_phage_PAS50';'Proteus_unclassified';'Providencia_unclassified';'Pseudomonas_phage_F116';'Pseudomonas_phage_MP38';'Pseudomonas_unclassified';'Ralstonia_unclassified';'Roseburia_unclassified';'Roseolovirus_unclassified';'Rothia_unclassified';'Saccharomyces_cerevisiae';'Salmonella_phage_FSL_SP_004';'Scardovia_unclassified';'Slackia_unclassified';'Staphylococcus_phage_PVL';'Streptococcus_mitis_oralis_pneumoniae';'Streptococcus_phage_Dp_1';'Streptococcus_phage_EJ_1';'Subdoligranulum_unclassified';'Thermus_unclassified';'Thiomonas_unclassified';'Turicibacter_unclassified';'Veillonella_unclassified';'Weissella_unclassified';'Yersinia_phage_L_413C';'Zika_virus'};
[considered,~] = setdiff(abundance(2:end,1),notConsidered);

% get percentage of captured taxa
sharedSpecies1=intersect(agoraSpecies1,considered);
sharedSpecies2=intersect(agoraSpecies2,considered);
agora1Captured=length(sharedSpecies1)/length(considered);
agora2Captured=length(sharedSpecies2)/length(considered);

calculatedCoverage{2,1}='CRC_Yachida2019';
calculatedCoverage{2,2}=agora1Captured;
calculatedCoverage{2,3}=agora2Captured;

% get percentage of captured abundance after cutting out unclassified
[C,IA]=intersect(abundance(:,1),notConsidered);
abundanceReduced=abundance;
abundanceReduced(IA,:)=[];

% get taxa covered by AGORA 1.03
[C,IA]=setdiff(abundanceReduced(:,1),agoraSpecies1,'stable');
abundanceAGORA1=abundanceReduced;
abundanceAGORA1(IA(2:end),:)=[];

% get taxa covered by AGORA 2
[C,IA]=setdiff(abundanceReduced(:,1),agoraSpecies2,'stable');
abundanceAGORA2=abundanceReduced;
abundanceAGORA2(IA(2:end),:)=[];

for i=2:size(abundance,2)
    data(i-1,1)=sum(cell2mat(abundance(2:end,i)));
    data(i-1,2)=sum(cell2mat(abundanceReduced(2:end,i)));
    data(i-1,3)=sum(cell2mat(abundanceAGORA1(2:end,i)));
    data(i-1,4)=sum(cell2mat(abundanceAGORA2(2:end,i)));
end
for i=1:length(data)
    reducedCaptured(i,1)=data(i,2)/data(i,1);
    agora1Captured(i,1)=data(i,3)/data(i,2);
    agora2Captured(i,1)=data(i,4)/data(i,2);
end
captureAfterReduction=[num2str(mean(reducedCaptured)) ' +/- ' num2str(std(reducedCaptured))];
calculatedCoverage{2,4}=[num2str(mean(agora1Captured)) ' +/- ' num2str(std(agora1Captured))];
calculatedCoverage{2,5}=[num2str(mean(agora2Captured)) ' +/- ' num2str(std(agora2Captured))];

%% calculate coverage of COMBO/PLEASE data

ibdFolder = [filesep 'Users' filesep 'mspg' filesep 'National University of Ireland, Galway/Group_MSP - Documents/Microbiome' filesep 'IBD_COMBO_PLEASE'];
pleaseData = readtable([ibdFolder filesep 'PLEASE' filesep 'S_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result.xlsx'], 'ReadVariableNames', false);
pleaseData = table2cell(pleaseData);
comboData = readtable([ibdFolder filesep 'COMBO' filesep 'S_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result.xlsx'], 'ReadVariableNames', false);
comboData = table2cell(comboData);

combinedData = pleaseData;
colSize=size(combinedData,2)+1;
combinedData(1,colSize:colSize+size(comboData,2)-2)=comboData(1,2:end);
combinedData(2:end,colSize:colSize+size(comboData,2)-2)={'0'};
for i=2:size(comboData,1)
    findRow=find(strcmp(combinedData(:,1),comboData{i,1}));
    if ~isempty(findRow)
        combinedData(findRow,colSize:colSize+size(comboData,2)-2)=comboData(i,2:end);
    else
        newRow=size(combinedData,1)+1;
        combinedData{newRow,1}=comboData{i,1};
        combinedData(newRow,2:end)={'0'};
        combinedData(newRow,colSize:colSize+size(comboData,2)-2)=comboData(i,2:end);
    end
end
rename={'Actinomyces_cardiffensis','Schaalia_cardiffensis';'Actinomyces_meyeri','Schaalia_meyeri';'Actinomyces_odontolyticus','Schaalia_odontolytica';'Actinomyces_turicensis','Schaalia_turicensis';'Clostridium_nexile','Tyzzerella_nexilis';'Eubacterium_hallii','Anaerobutyricum_hallii';'Ruminococcus_gnavus','Blautia_gnavus';'Ruminococcus_torques','Blautia_torques';'Propionibacterium_acnes','Cutibacterium_acnes';'Propionibacterium_avidum','Cutibacterium_avidum';'Acinetobacter_baumanni','Acinetobacter_baumannii';'Clostridium_difficile','Clostridioides_difficile';'Enterobacter_aerogenes','Klebsiella_aerogenes';'Streptococcus_bovis','Streptococcus_equinus';'Clostridium_bartlettii','Intestinibacter_bartlettii';'Clostridium_hathewayi','Hungatella_hathewayi';'Gemella_moribillum','Gemella_morbillorum';'Ruminococcus_obeum','Blautia_obeum';'Prevotella_tannerae','Alloprevotella_tannerae';'Eubacterium_saburreum','Lachnoanaerobaculum_saburreum'};
for i=1:length(rename)
    combinedData(:,1)=strrep(combinedData(:,1),rename{i,1},rename{i,2});
end
% needs to be saved manually as csv file-too big to export

% create the mapping

% AGORA 1.03
fileDir = fileparts(which('ReactionTranslationTable.txt'));
cd(fileDir);
[translatedAbundances1,normalizedAbundances1,unmappedRows1]=translateMetagenome2AGORA([ibdFolder filesep 'COMBO_PLEASE_data.csv'],'Species');
% AGORA 2.0
[translatedAbundances2,normalizedAbundances2,unmappedRows2]=translateMetagenome2AGORA2([ibdFolder filesep 'COMBO_PLEASE_data.csv'],'Species');

% get percentage of captured taxa

agora1Captured=(size(translatedAbundances1,1)-1)/(size(combinedData,1)-1);
agora2Captured=(size(translatedAbundances2,1)-1)/(size(combinedData,1)-1);

calculatedCoverage{3,1}='IBD_Lewis2015';
calculatedCoverage{3,2}=agora1Captured;
calculatedCoverage{3,3}=agora2Captured;

for i=2:size(combinedData,2)
    for j=2:size(combinedData,1)
        if ~isnumeric(combinedData{j,i})
            combinedData{j,i}=str2double(combinedData{j,i});
        end
    end
    data(i-1,1)=sum(cell2mat(combinedData(2:end,i)));
    for j=2:size(translatedAbundances1,1)
        if ~isnumeric(translatedAbundances1{j,i})
            translatedAbundances1{j,i}=str2double(translatedAbundances1{j,i});
        end
    end
    data(i-1,2)=sum(cell2mat(translatedAbundances1(2:end,i)));
    for j=2:size(translatedAbundances2,1)
        if ~isnumeric(translatedAbundances2{j,i})
            translatedAbundances2{j,i}=str2double(translatedAbundances2{j,i});
        end
    end
    data(i-1,3)=sum(cell2mat(translatedAbundances2(2:end,i)));
end
for i=1:length(data)
    agora1Captured(i,1)=data(i,2)/data(i,1);
    agora2Captured(i,1)=data(i,3)/data(i,1);
end
calculatedCoverage{3,4}=[num2str(mean(agora1Captured)) ' +/- ' num2str(std(agora1Captured))];
calculatedCoverage{3,5}=[num2str(mean(agora2Captured)) ' +/- ' num2str(std(agora2Captured))];
