function summaryMerged = collectPhylumStats(microbiotaPath, taxonomyPath, microbiotaWbmPath)
% How did mapping change the phylum level relative abundances?

% Raw microbiota relative abundances with processed names
microbiota = readtable(microbiotaPath);
microbiota = renamevars(microbiota,'Taxon','Species');

% Species-phylum mapping after species renaming
% Load Taxonomy information
opts = detectImportOptions(taxonomyPath);
opts.SelectedVariableNames = "Taxon";
taxonomyInfo = readmatrix(taxonomyPath, opts);
taxonomyInfo = append(taxonomyInfo,'|');

% Information on which microbial species were mapped on AGOR2+APOLLO

% Split taxonomic data into multiple columns
levels = {'Kingdom','Phylum','Clade','Order','Family','Genus','Species'};
levelAbbreviations = ['k','p','c','o','f','g','s'];
expresFun = @(x) ['(?<=',x,'__)(.*?)(?=\|)']; % Get the taxonomy information of interest
regexFun = @(x) string(regexp(taxonomyInfo,expresFun(x),'tokens','once')); % Extract matches in regex
taxaInfo = arrayfun(regexFun, levelAbbreviations,'UniformOutput',false); % Run regex for all levels
taxaInfo = array2table(cellstr(horzcat(taxaInfo{:})),'VariableNames',levels); % Generate table

% Select phylum infotaxaInfo
taxaInfo = taxaInfo(:,{'Species','Phylum'});

% Map species names on taxaInfo. 
% microbiota = microbiota(:,{'Species'});

mergedTable = outerjoin(microbiota,taxaInfo,"Keys","Species","Type","left",'MergeKeys',true);
mergedTable = movevars(mergedTable,'Phylum','After','Species');

% 24 microbial species could not be matches with their phylum information.
% This information will now be added manually using the NCBI taxonomy
% browser (accessed on April 7, 2025):
% https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

% Apply changes
mergedTable.Phylum(contains(mergedTable.Species,{'Actinomyces','Actinobaculum'})) = {'Actinobacteria'};
mergedTable.Phylum(contains(mergedTable.Species,{'Prevotella','Bacteroides','Dysgonomonas'})) = {'Bacteroidetes'};
mergedTable.Phylum(contains(mergedTable.Species,{'Blautia','Clostridiales','Enterococcus','Eubacterium','Lachnospiraceae','Megasphaera','Oscillibacter','Parvimonas','Stoquefichus'})) = {'Firmicutes'};
mergedTable.Phylum(contains(mergedTable.Species,{'Methanomassiliicoccales'})) = {'Euryarchaeota'};

% Identify which microbial species in the raw data were mapped.
% Load processed and mapped species relative abundances 
microbiotaWBM = readtable(microbiotaWbmPath,'ReadRowNames',true);
microbiotaWBM.SumOfTaxa = [];

% Process microbiota data
microbiotaWBM = renamevars(rows2vars(microbiotaWBM),'OriginalVariableNames','Species');
microbiotaWBM.Species = replace(microbiotaWBM.Species,'_',' ');

% Repair names in the mapped taxa to ensure overlap
mergedTable.Species(matches(mergedTable.Species,'Ruminococcus champanellensis')) = {'Ruminococcus chamellensis'};
mergedTable.Species(matches(mergedTable.Species,'Companilactobacillus farciminis')) = {'Comilactobacillus farciminis'};

% Add information on which microbial species were mapped onto the WBMs
mergedTable.mapped = matches(mergedTable.Species,microbiotaWBM.Species);
mergedTable = movevars(mergedTable,'mapped','After','Phylum');
%%
sumStats = cell(1,3);
sumStats{1} = getPhylumSummary(mergedTable); % Summary statistics on all microbes
sumStats{2} = getPhylumSummary(mergedTable(mergedTable.mapped==1,:)); % Summary statistics on mapped microbes
sumStats{3} = getPhylumSummary(mergedTable(mergedTable.mapped==0,:)); % Summary statistics on unmapped microbes

% Combine tables 
summaryMerged = outerjoin(sumStats{1},sumStats{2},"Keys","Phylum","Type","left",'MergeKeys',true); % Combine total with mapped
summaryMerged = outerjoin(summaryMerged,sumStats{3},"Keys","Phylum","Type","left",'MergeKeys',true); % Combine total+mapped with unmapped

% Rename variables
newNames = {'Phylum','Mean_premapped','SD_premapped','Mean_mapped','SD_mapped','Mean_unmapped','SD_unmapped'};
summaryMerged = renamevars(summaryMerged,summaryMerged.Properties.VariableNames,newNames);
end


function sumStats = getPhylumSummary(input)
% Calculate the phylum-level read counts to find the raw phylum-level
% relative abundances.
phylumTable = removevars(input,{'Species','mapped'});

% Calculate the phylum level relative abundances
relAbun = phylumTable{:,2:end};
phyla = phylumTable.Phylum;

% Find the phylum groups rows
[rows,groups] = findgroups(phyla);

% Calculate group-wise summed relative abundances
pRA = arrayfun(@(x) splitapply(@sum,relAbun(:,x),rows),1:width(relAbun),'UniformOutput',false);
pRA = horzcat(pRA{:}); % Create matrix
pRA = pRA ./ sum(pRA); % Renormalise
sumStats = table(groups, mean(pRA,2) ,std(pRA,[],2),'VariableNames',{'Phylum','Mean','SD'});
end

