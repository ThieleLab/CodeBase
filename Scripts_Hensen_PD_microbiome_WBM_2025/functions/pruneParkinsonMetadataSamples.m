function [prunedMetada, sampleNumbers] = pruneParkinsonMetadataSamples(metadata)
% This function defines the variables in the metadata of interest and
% removes all samples that miss one or more samples in the metadata. Sample
% summary statistics are also generated.
%
% USAGE:
%   metadata = pruneParkinsonMetadataSamples(metadataPath)
%
% INPUT:
%   metadataPath: string
%       Path to the metadata table in a readable format (e.g., .xlsx or .csv).
%
% OUTPUT:
%   metadataPathRmSamplesPath:  Path to pruned metadata
%   metadataRmSummaryPath:      Path to summary table with the number of
%                               removed samples.
%
% AUTHORS:
%   Tim Hensen, December 2024

% Create statistics table with sample numbers, before and after sample
% pruning.

sampleNumbers = table('Size',[3,4],...
    'VariableTypes',{'double','double','double','double'},...
    'VariableNames',{'Total','PD','Control','Removed'},...
    'RowNames',	{'Original','Removed swabs','Removed missing data'});

% Populate the row with original sample numbers
sampleNumbers{"Original","Total"} = size(metadata,1);
sampleNumbers{"Original","PD"} = sum(matches(metadata.Case_status,'PD'));
sampleNumbers{"Original","Control"} = sum(matches(metadata.Case_status,'Control'));


% Step 1: Find and remove all individuals with swab method collection
metadataOri = metadata;
metadata(matches(metadata.collection_method,"swab"),:) = [];

% Populate the row with new sample numbers
sampleNumbers{"Removed swabs","Total"} = size(metadata,1);
sampleNumbers{"Removed swabs","PD"} = sum(matches(metadata.Case_status,'PD'));
sampleNumbers{"Removed swabs","Control"} = sum(matches(metadata.Case_status,'Control'));
sampleNumbers{"Removed swabs","Removed"} = size(metadataOri,1) - size(metadata,1);


% Step 2: Define all variables not in the confounder analysis

% PD/NHC, total sequence count, intake (yes/no) of alcohol, laxatives, 
% probiotics, antihistamines, depression-anxiety-mood medication, pain medication, and sleep-aid. 
variablesToCheck = {'ID','Case_status','Sex','Age_at_collection',...
    'Do_you_drink_alcohol','Laxatives','Probiotic','Antihistamines',...
    'Depression_anxiety_mood_med','Pain_med','Sleep_aid'};


% Step 3: Find and remove all samples with missing data in one or more
% confounders. 
metadataOri = metadata;
for i = 1:length(variablesToCheck)
    variable = metadata.(variablesToCheck{i});
    if isnumeric(variable)
        noInfo = isnan(variable);
    else
        noInfo = cellfun(@isempty,variable);
    end
    metadata(noInfo,:) = [];
end

% Define function output
prunedMetada = metadata;

% Populate the row with final sample numbers
sampleNumbers{"Removed missing data","Total"} = size(metadata,1);
sampleNumbers{"Removed missing data","PD"} = sum(matches(metadata.Case_status,'PD'));
sampleNumbers{"Removed missing data","Control"} = sum(matches(metadata.Case_status,'Control'));
sampleNumbers{"Removed missing data","Removed"} = size(metadataOri,1) - size(metadata,1);

end