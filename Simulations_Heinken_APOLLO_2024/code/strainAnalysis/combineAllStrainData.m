% combine model properties for Pasolli and Almeida reconstructions into one
% file

mkdir([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'CombinedProperties' filesep 'ReactionMetabolitePresence'])
mkdir([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'CombinedProperties' filesep 'ComputedFluxes'])

datasets={
    ['ReactionMetabolitePresence' filesep 'ReactionPresence_']
    ['ReactionMetabolitePresence' filesep 'MetabolitePresence_']
    ['ComputedFluxes' filesep 'InternalProduction_']
    ['ComputedFluxes' filesep 'UptakeSecretion_']
    };

for i=1:length(datasets)
    feats_combined={};

    fid = fopen([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Almeida' filesep datasets{i} 'Almeida_refined.txt']);
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

    fid = fopen([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Pasolli' filesep datasets{i} 'Pasolli_refined.txt']);
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
    fid = fopen([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Almeida' filesep datasets{i} 'Almeida_refined.txt']);
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
    save([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'CombinedProperties' filesep datasets{i} 'combined.mat'],'data_combined','-v7.3');

    fid = fopen([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'refinedModelProperties_Pasolli' filesep datasets{i} 'Pasolli_refined.txt']);
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

    save([rootDir filesep 'data' filesep 'analysis_ModelProperties' filesep 'CombinedProperties' filesep datasets{i} 'combined.mat'],'data_combined','-v7.3');
end
