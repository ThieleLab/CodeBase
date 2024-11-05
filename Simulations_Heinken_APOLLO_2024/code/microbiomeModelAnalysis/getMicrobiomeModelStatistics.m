
% get the statistics of all microbiome models

clear all
rootDir = pwd;

mkdir([rootDir filesep 'results' filesep 'microbiomes'])

statistics = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'ModelStatistics.csv']);
statistics{1,4} = 'Microbes';

normCoverage = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_MicrobiomeModels' filesep 'Combined_data' filesep 'normalizedAbundance.csv']);

for i=2:size(normCoverage,2)
    findMod = find(strcmp(statistics(:,1),normCoverage{1,i}));
    abun = cell2mat(normCoverage(2:end,i));
    statistics{findMod,4} = length(find(abun>0.00001));
end

cell2csv([rootDir filesep 'results' filesep 'microbiomes' filesep 'MicrobiomeModelStatistics.csv'],statistics)