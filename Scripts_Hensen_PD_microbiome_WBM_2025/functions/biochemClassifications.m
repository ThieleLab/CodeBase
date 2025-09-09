function [biomchemClass,capturedSubsystems] = biochemClassifications(metaboliteClassPath)
% DESCRIPTION:
%   This function processes metabolite classification data to generate a summary of the distribution
%   of metabolites across biochemical categories (e.g., superclasses, classes, subclasses, Level 5 classes)
%   and metabolic subsystems. Additionally, it computes the fraction of metabolic subsystems captured
%   by the given set of metabolites relative to the VMH database.
%
% USAGE:
%   [biomchemClass, capturedSubsystems] = biochemClassifications(metaboliteClassPath)
%
% INPUTS:
%   metaboliteClassPath: String specifying the path to the input file containing
%                        metabolite classifications. The file should be in a tabular format
%                        with at least the following columns:
%                        - 'superclass': The superclass classification of each metabolite.
%                        - 'class': The class-level classification of each metabolite.
%                        - 'subclass': The subclass-level classification of each metabolite.
%                        - 'Level 5': Optional column representing additional classification levels.
%                        - 'Subsystem': Associated metabolic subsystem for each metabolite.
%   saveDir: string
%       Path to the directory where output files (table and Venn diagram) will be saved.
%
% OUTPUTS:
%   biomchemClass: A structure containing the biochemical classifications of metabolites:
%                  - biomchemClass.superClasses: Table summarizing the count of metabolites in each superclass.
%                  - biomchemClass.classes: Table summarizing the count of metabolites in each class.
%                  - biomchemClass.subclasses: Table summarizing the count of metabolites in each subclass.
%                  - biomchemClass.level5classes: Table summarizing the count of metabolites in each Level 5 classification.
%                  - biomchemClass.subsystems: Table summarizing the count of metabolites in each subsystem.
%
%   capturedSubsystems: Fraction of unique metabolic subsystems covered by the input data compared to
%                       the subsystems in the Virtual Metabolic Human (VMH) database. This is calculated
%                       as the ratio of subsystems present in the input data to the total number of unique
%                       subsystems in the VMH database.
%
% AUTHOR:
%   Tim Hensen, January 2025

% Load metabolite annotations
classifications = readtable(metaboliteClassPath,'VariableNamingRule','preserve');

% Find the number of groups and the number of metabolites in each group
biomchemClass = struct;
biomchemClass.classifications = classifications;

% Superclasses
[groups,names] = findgroups(classifications.superclass);
biomchemClass.superClasses = table(names,histcounts(groups)','VariableNames',{'Super class','Metabolites'});

% Classes
[groups,names] = findgroups(classifications.class);
biomchemClass.classes = table(names,histcounts(groups)','VariableNames',{'Class','Metabolites'});

% Sub classes
[groups,names] = findgroups(classifications.subclass);
biomchemClass.subclasses = table(names,histcounts(groups)','VariableNames',{'Sub class','Metabolites'});

% Level 5 classification
[groups,names] = findgroups(classifications.("Level 5"));
biomchemClass.level5classes = table(names,histcounts(groups)','VariableNames',{'Level 5 class','Metabolites'});

% Associated subsystems
[groups,names] = findgroups(classifications.Subsystem);
biomchemClass.subsystems = table(names,histcounts(groups)','VariableNames',{'Subsystem','Metabolites'});

% Calculate the fraction of total metabolic subsystems that are captured by
% the 116 metabolites.
database = loadVMHDatabase;
metabolites = database.metabolites;
metabolites(cellfun(@isempty,metabolites(:,13)),:)=[];

capturedSubsystems = height(biomchemClass.subsystems) /  numel(unique(metabolites(:,13)));

% Save biochemical classifications to results
% Write each field to a separate sheet
% filename = [saveDir filesep 'metaboliteBiochemAnnotations.xlsx'];
% fieldNames = string(fieldnames(biomchemClass));
% arrayfun(@(x) writetable(biomchemClass.(x), filename, 'Sheet', x), fieldNames, 'UniformOutput', false);

end

