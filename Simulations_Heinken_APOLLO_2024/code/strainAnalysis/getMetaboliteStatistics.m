
% get statistics on the number of produced and secreted metabolites

resources={
    'Pasolli'
    'Almeida'
    };

stats = {'Unique metabolites','Secreted metabolites','Consumed metabolites','Internally produced metabolites'};
for i=1:length(stats)
    stats{2,i}={};
end

for i=1:size(resources,1)
    
    % secreted metabolites
    data = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_' resources{i,1} filesep 'ComputedFluxes' filesep 'secretionFluxes_' resources{i,1} '_refined.txt']);
    
    for j=2:size(data,2)
        if any(cell2mat(data(2:end,j))>0.000001)
            stats{2,2}=union(stats{2,2},data{1,j});
        end
    end

    % consumed metabolites
    data = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_' resources{i,1} filesep 'ComputedFluxes' filesep 'uptakeFluxes_' resources{i,1} '_refined.txt']);
    for j=2:size(data,2)
        if any(cell2mat(data(2:end,j))<-0.000001)
            stats{2,3}=union(stats{2,3},data{1,j});
        end
    end

     % internally produced metabolites
    data = readInputTableForPipeline([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_' resources{i,1} filesep 'ComputedFluxes' filesep 'InternalProduction_' resources{i,1} '_refined.txt']);
    for j=2:size(data,2)
        stats{2,1}=union(stats{2,1},data{1,j});
        try
            if any(cell2mat(data(2:end,j))>0.000001) 
                stats{2,4}=union(stats{2,4},data{1,j});
            end
        end
    end
    save([rootDir filesep 'results' filesep 'strains' filesep 'MetaboliteStats.mat'],'stats')
end
