function mergedTaxa = collectAbundanceStats(paths, fileNames)

% Load the pre-mapped, mapped present, and mapped absent abundance data
abundances = cellfun(@(x) readtable(fullfile(paths.species,x)),fileNames,'UniformOutput',false);

% Remove pan indications from mapped microbial species
abundances{2}.Taxon = erase(abundances{2}.Taxon,'pan');

% Repair incorrect names in the mapped taxa to ensure overlap
oldNames = {'Ruminococcus chamellensis','Comilactobacillus farciminis'};
nameRepairs = {'Ruminococcus champanellensis','Companilactobacillus farciminis'}';
abundances{2}.Taxon(matches(abundances{2}.Taxon, oldNames)) = nameRepairs; 

% Remove columns
rmCols = @(x) removevars(x,{'maximum','minimum','non_zero_count'});
abundances = cellfun(rmCols,abundances,'UniformOutput',false);

% Annotate columns for future concatenation
renameCols = @(x,y) renamevars(x,2:width(x), y);
newNames = {{'mean_premapped','SD_premapped'},{'mean_mapped','SD_mapped'},{'mean_unmapped','SD_unmapped'}};
abundances = cellfun(renameCols,abundances,newNames,'UniformOutput',false);

% Add labels on which microbial species are mapped
abundances{2}.mapped = true(height(abundances{2}),1);

% Combine tables 
mergedTaxa = outerjoin(abundances{1},abundances{3},'Keys','Taxon','MergeKeys',true);
mergedTaxa = outerjoin(mergedTaxa,abundances{2},'Keys','Taxon','MergeKeys',true);

% Rename columns
mergedTaxa = movevars(mergedTaxa,'mapped','After','Taxon');
mergedTaxa = renamevars(mergedTaxa,'Taxon','Microbial species');
mergedTaxa = renamevars(mergedTaxa,'mapped','Mapped on AGORA2+APOLLO');
end