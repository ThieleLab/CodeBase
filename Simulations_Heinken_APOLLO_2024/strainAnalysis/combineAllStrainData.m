% combine strain-level data for 150k and 90k

mkdir(['Model_properties_analysis' filesep 'CombinedProperties' filesep 'ReactionMetabolitePresence'])
mkdir(['Model_properties_analysis' filesep 'CombinedProperties' filesep 'ComputedFluxes'])

datasets={
    ['ReactionMetabolitePresence' filesep 'ReactionPresence_']
    ['ReactionMetabolitePresence' filesep 'MetabolitePresence_']
    ['ComputedFluxes' filesep 'InternalProduction_']
    ['ComputedFluxes' filesep 'UptakeSecretion_']
    };

for i=1:length(datasets)
    feats_combined={};

    fid = fopen(['Model_properties_analysis' filesep '90k_properties' filesep datasets{i} '90k_refined.txt']);
    tline = fgetl(fid);
    tlines = cell(0,1);
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
    headers=strsplit(tlines{1,1},'	');
    clear tlines

    feats_combined=union(feats_combined,headers(2:end));

    fid = fopen(['Model_properties_analysis' filesep '150k_properties' filesep datasets{i} '150k_refined.txt']);
    tline = fgetl(fid);
    tlines = cell(0,1);
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
    headers=strsplit(tlines{1,1},'	');
    clear tlines
    feats_combined=union(feats_combined,headers(2:end));

    data_combined = {''};
    for j=1:length(feats_combined)
        data_combined{1,1} = [data_combined{1,1} '_spl_' feats_combined{j}];
    end
    % find indices in original tables for both datasets
    fid = fopen(['Model_properties_analysis' filesep '90k_properties' filesep datasets{i} '90k_refined.txt']);
    tline = fgetl(fid);
    tlines = cell(0,1);
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
    headers = strsplit(tlines{1,1},'	');
    combined_headers = strsplit(data_combined{1,1},'_spl_');
    [~,origFeatInd,combFeatInd] = intersect(headers,combined_headers,'stable');

    cnt=2;
    for j=2:size(tlines,1)
        j
        % retrieve the data
        data_tmp = strsplit(tlines{j,1},'	');
        data_combined{cnt,1}=data_tmp{1,1};
        for k=2:length(origFeatInd)
            data_combined{cnt,1}=[data_combined{cnt,1} '_spl_' data_tmp{1,origFeatInd(k,1)}];
        end
        cnt=cnt+1;
    end
    save(['Model_properties_analysis' filesep 'CombinedProperties' filesep datasets{i} 'combined.mat'],'data_combined','-v7.3');

    fid = fopen(['Model_properties_analysis' filesep '150k_properties' filesep datasets{i} '150k_refined.txt']);
    tline = fgetl(fid);
    tlines = cell(0,1);
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
    headers=strsplit(tlines{1,1},'	');
    combined_headers = strsplit(data_combined{1,1},'_spl_');
    [~,origFeatInd,combFeatInd] = intersect(headers,combined_headers,'stable');

    for j=2:size(tlines,1)
        j
        % retrieve the data
        data_tmp = strsplit(tlines{j,1},'	');
        data_combined{cnt,1}=data_tmp{1,1};
        for k=2:length(origFeatInd)
            data_combined{cnt,1}=[data_combined{cnt,1} '_spl_' data_tmp{1,origFeatInd(k,1)}];
        end
        cnt=cnt+1;
    end

    save(['Model_properties_analysis' filesep 'CombinedProperties' filesep datasets{i} 'combined.mat'],'data_combined','-v7.3');
end

% export the combined data
for i=1:length(datasets)
    load(['Model_properties_analysis' filesep 'CombinedProperties' filesep datasets{i} 'combined.mat'])

    data = {};
    for j=1:size(data_combined,1)
        j
        data_tmp = strsplit(data_combined{j,1},'_spl_');
        for k=1:length(data_tmp)
            data{j,k}=data_tmp{k};
        end
    end
    writetable(cell2table(data),['Model_properties_analysis' filesep 'CombinedProperties' filesep datasets{i} 'combined'],'WriteVariableNames',false)
end
