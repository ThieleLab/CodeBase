function suplMetaboliteTable = collectMetInfo(paths_metWbmLocations,paths_metOntology)
% Combine metabolite WBM location information with chemical ontology
% information table.

% Load metabolite WBM info
metWbmLocations = readtable(paths_metWbmLocations,'VariableNamingRule','preserve');

% Rename variables
currName = {'Metabolite','blood','Diet','Microbiome'};
newName = {'VMH ID','Present in WBM blood compartment','Present in given diet','Microbial metabolite'};
metWbmLocations = renamevars(metWbmLocations,currName,newName);

% Filter on columns of interest
metWbmLocations = metWbmLocations(:,newName);

% Remove metabolites that were not analysed
metWbmLocations(metWbmLocations.("Present in WBM blood compartment")==0,:)=[];

% Read and process metabolite taxonomy data
metOntology = readtable(paths_metOntology,'VariableNamingRule','preserve');
metOntology = renamevars(metOntology,{'Metabolite','Common name','Subsystem'},{'VMH ID','Metabolite name','Metabolic subsystem'});

% Combine datasets
suplMetaboliteTable = outerjoin(metOntology,metWbmLocations,'Type','full','Keys','VMH ID','MergeKeys',true);
end