
% get statistics on the number of produced and secreted metabolites

versions = {
    '150k','150k_properties'
    '90k','90k_properties'
    };

stats = {'Unique metabolites','Secreted metabolites','Consumed metabolites','Internally produced metabolites'};
for i=1:length(stats)
    stats{2,i}={};
end

for i=1:size(versions,1)
    
    % secreted metabolites
    data = readInputTableForPipeline(['Model_properties_analysis' filesep versions{i,2} filesep 'Refined' filesep 'ComputedFluxes' filesep 'secretionFluxes_' versions{i,1} '_refined.txt']);
    
    for j=2:size(data,2)
        if any(cell2mat(data(2:end,j))>0.000001)
            stats{2,2}=union(stats{2,2},data{1,j});
        end
    end
    save(['Overview_Model_stats' filesep 'MetaboliteStats.mat'],'stats')

    % consumed metabolites
    data = readInputTableForPipeline(['Model_properties_analysis' filesep versions{i,2} filesep 'Refined' filesep 'ComputedFluxes' filesep 'uptakeFluxes_' versions{i,1} '_refined.txt']);
    for j=2:size(data,2)
        if any(cell2mat(data(2:end,j))<-0.000001)
            stats{2,3}=union(stats{2,3},data{1,j});
        end
    end
    save(['Overview_Model_stats' filesep 'MetaboliteStats.mat'],'stats')

     % internally produced metabolites
    data = readInputTableForPipeline(['Model_properties_analysis' filesep versions{i,2} filesep 'Refined' filesep 'ComputedFluxes' filesep 'InternalProduction_' versions{i,1} '_refined.txt']);
    for j=2:size(data,2)
        stats{2,1}=union(stats{2,1},data{1,j});
        try
            if any(cell2mat(data(2:end,j))>0.000001) 
                stats{2,4}=union(stats{2,4},data{1,j});
            end
        end
    end
    save(['Overview_Model_stats' filesep 'MetaboliteStats.mat'],'stats')
end
