
% get the statistics of all microbiome models

statistics = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Combined_data' filesep 'ModelStatisticsCombined.csv']);
statistics{1,4} = 'Microbes';

normCoverage = readInputTableForPipeline(['Analysis_microbiome_models' filesep 'Combined_data' filesep 'normalizedAbundanceCombined.csv']);

for i=2:size(normCoverage,2)
    findMod = find(strcmp(statistics(:,1),normCoverage{1,i}));
    abun = cell2mat(normCoverage(2:end,i));
    statistics{findMod,4} = length(find(abun>0.00001));
end

cell2csv(['Overview_Model_stats' filesep 'MicrobiomeModelStatistics.csv'],statistics)