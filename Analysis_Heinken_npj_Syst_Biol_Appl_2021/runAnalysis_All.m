%%
% This script reproduces the analysis performed in Heinken et al, 
% "Metabolic modelling reveals broad changes in gut microbial metabolism in 
% inflammatory bowel disease patients with dysbiosis."
% The script relies on functions in the COBRA Toolbox
% (https://github.com/opencobra/cobratoolbox).
%   - Almut Heinken, 03/2021

%% Define necessary inputs
% download and set path to AGORA reconstructions
system('curl -LJO https://github.com/VirtualMetabolicHuman/AGORA/archive/master.zip')
unzip('AGORA-master')
agoraPath = [pwd filesep 'AGORA-master' filesep 'CurrentVersion' filesep 'AGORA_1_03' filesep' 'AGORA_1_03_mat'];

% download and extract the necessary input file
websave('inputFiles','https://7d0178ba-85cc-458f-bb4a-8b987a027734.filesusr.com/archives/f63236_59cc1a0a3454476d873c39ae86418749.zip')
unzip('inputFiles')

% download and extract the models and fluxes generated in the study
websave('Computed_fluxes','https://7d0178ba-85cc-458f-bb4a-8b987a027734.filesusr.com/archives/f63236_9bb5aae408074167b9335ae64f5231cb.zip')
unzip('Computed_fluxes')

mkdir('Constrained_models')
cd('Constrained_models')
websave('Constrained_models_1','https://7d0178ba-85cc-458f-bb4a-8b987a027734.filesusr.com/archives/f63236_6beb8e3e2cf14706ad10ed275429e720.zip')
unzip('Constrained_models_1')
websave('Constrained_models_2','https://7d0178ba-85cc-458f-bb4a-8b987a027734.filesusr.com/archives/f63236_a23b350248dc487b9e83985730ce1e8e.zip')
unzip('Constrained_models_2')
cd ..

% path to the file with normalized abundances for COMBO/PLEASE samples
abunFilePath=[pwd filesep 'inputFiles' filesep 'normalizedCoverage_IBD.csv'];

% path to file where net uptake and secretion fluxes will be saved
netFluxesPath = [pwd filesep 'NetFluxes'];
mkdir(netFluxesPath)

% path to where strain to metabolite contributions are saved
contrPath = [pwd filesep 'MicrobeContributions'];
mkdir(contrPath)

% parth to where strain to metabolite contribution profiles for each
% metabolite are saved
profilePath = [pwd filesep 'MetaboliteProfiles'];
mkdir(profilePath)

% path to where correlations between net production fluxes and taxon
% abundances are saved
corrPath = [pwd filesep 'Correlations'];
mkdir(corrPath)

% path to where reaction and subsystem abundances are saved
rxnsPath = [pwd filesep 'ReactionAbundances'];
mkdir(rxnsPath)

% path to file with sample information
metadataPath = [pwd filesep 'inputFiles' filesep 'metadata_IBD.csv'];

% load taxonomical information on AGORA organisms
infoFile = readtable('AGORA_infoFile.xlsx', 'ReadVariableNames', false);
infoFile = table2cell(infoFile);

% load table with metabolite information
metaboliteInfo=table2cell(readtable([pwd filesep 'inputFiles' filesep 'MetaboliteInformation.csv'], 'ReadVariableNames', false));

% number of cores dedicated for parallelization
numWorkers = 4;

addpath('scripts')

%% run all scripts in order

%% Create microbiome models with dietary constraints
% NOTE: due to ongoing changes in the COBRA Toolbox, newly created models
% may be slightly different from the ones used for the original analysis.
% We provide the script to build the models, but to reproduce the results,
% we provide the original microbiome models as a download.

% To build personalized microbiome models:
% runMgPipe

% To download the originally used microbiome models:
modelPath=[pwd filesep 'Constrained_models'];

%% Compute metabolite secretion profiles and strain to metabolite contributions
% This step is done in Julia through the COBRA.jl package.
% See the file startDistributedFBA.txt for instructions how to run the 
% computation in COBRA.jl. Note that the customized Julia script may not be
% compatible with newer versions of Julia/COBRA.jl.

%% Extraction of computed fluxes
addpath([pwd filesep 'scripts']);

%% Net uptake and secretion
fluxPath=[pwd filesep 'Computed_fluxes' filesep 'Net_uptake_secretion_fluxes'];
extractNetFluxes

%% create violin plots for net secretion fluxes
metadata = readtable(metadataPath, 'ReadVariableNames', false);
metadata = table2cell(metadata);

violinPath = [pwd filesep 'ViolinPlots'];
mkdir(violinPath)

sampleData = readtable([netFluxesPath filesep 'netProductionFluxes_Exported.csv'], 'ReadVariableNames', false);
sampleData = table2cell(sampleData);

cd(violinPath)
makeViolinPlots(sampleData, metadata, 'unit', 'mmol/person/day')
cd ..

%% Extract strain to metabolite contributions
contrFluxPath=[pwd filesep 'Computed_fluxes' filesep 'Strain_metabolite_contributions'];
extractMicrobeContributions

%% extract strain contribution profiles separately for each metabolite
extractMetaboliteProfiles

%% Calculate reaction presence, reaction abundance, and subsystem abundance
calculateReactionAbundances

%% Calculate correlations between taxon abundances and net fluxes
correlateFluxesWithTaxa

%% Perform statistical analysis of the results
statPath = [pwd filesep 'Statistics'];

files={
    [netFluxesPath filesep 'netProductionFluxes.csv'],'netProductionFluxes'
    [netFluxesPath filesep 'netUptakeFluxes.csv'],'netUptakeFluxes'
    [contrPath filesep 'MicrobeContributions_Fluxes.csv'],'MicrobeContributions_Fluxes'
    [rxnsPath filesep 'ReactionPresence.csv'],'ReactionPresence'
    [rxnsPath filesep 'ReactionAbundance.csv'],'ReactionAbundance'
    [rxnsPath filesep 'ReactionAbundance_Phylum.csv'],'ReactionAbundance_Phylum'
    [rxnsPath filesep 'ReactionAbundance_Genus.csv'],'ReactionAbundance_Genus'
    [rxnsPath filesep 'SubsystemAbundance.csv'],'SubsystemAbundance'
    };

% significant differences between IBD and healthy
metadataPath = [pwd filesep 'inputFiles' filesep 'metadata_IBD.csv'];
metadata = table2cell(readtable(metadataPath, 'ReadVariableNames', false));
metadata(:,2)=strrep(metadata(:,2),'IBD_dysbiotic','IBD');
metadata(:,2)=strrep(metadata(:,2),'IBD_nondysbiotic','IBD');

for i=1:size(files,1)
    sampleData = table2cell(readtable(files{i,1}, 'ReadVariableNames', false));
    Statistics = performStatisticalAnalysis(sampleData',metadata);
    writetable(cell2table(Statistics),[statPath filesep files{i,2} '_Statistics_Healthy_vs_IBD'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
end

% significant differences between dysbiotic and nondysbiotic IBD
metadataPath = [pwd filesep 'inputFiles' filesep 'metadata_IBD.csv'];
metadata = table2cell(readtable(metadataPath, 'ReadVariableNames', false));
metadata(find(strcmp(metadata(:,2),'Healthy')),:)=[];

for i=1:size(files,1)
    sampleData = table2cell(readtable(files{i,1}, 'ReadVariableNames', false));
    % remove healthy samples
    [C,IA]=setdiff(sampleData(1,:),metadata(:,1),'stable');
    sampleData(:,IA(2:end))=[];
    Statistics = performStatisticalAnalysis(sampleData',metadata);
    writetable(cell2table(Statistics),[statPath filesep files{i,2} '_Statistics_IBD_nondysbiotic_vs_dysbiotic'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
end

% get all features that were significantly different
for i=1:size(files,1)
    sampleData = table2cell(readtable(files{i,1}, 'ReadVariableNames', false));
    
    % get significant features
    statistics = table2cell(readtable([statPath filesep files{i,2} '_Statistics_Healthy_vs_IBD.txt'], 'ReadVariableNames', false));
    sigFeat=statistics(find(strcmp(statistics(:,4),'1')),1);
    statistics = table2cell(readtable([statPath filesep files{i,2} '_Statistics_IBD_nondysbiotic_vs_dysbiotic.txt'], 'ReadVariableNames', false));
    sigFeat=union(sigFeat,statistics(find(strcmp(statistics(:,4),'1')),1));
    
    % delete all features that were nonsignificant
    [C,IA]=setdiff(sampleData(:,1),sigFeat,'stable');
    sampleData(IA(2:end),:)=[];
    writetable(cell2table(sampleData),[statPath filesep files{i,2} '_SignificantFeatures'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
end
