% This script extracts the computed strain-level contributions to 
% metabolites for the personalized microbiome models generated for the 
% PLEASE/COMBO cohort.
% Almut Heinken, 03/21

metadata=table2cell(readtable(metadataPath, 'ReadVariableNames', false));

ContributionFluxes{1,2}='Strain';
ContributionFluxes{1,3}='Metabolite';

taxa=infoFile(1,3:8);
for i=1:length(taxa)
    ContributionFluxes{1,i+3}=taxa{i};
end

for i=2:size(metadata,1)
    ContributionFluxes{1,i+8}=metadata{i,1};
end

cnt=2;
for i=2:size(metadata,1)
    i
    % load model
    load([modelPath filesep metadata{i,1} '.mat']);
    % load corresponding contribution fluxes
    load([contrFluxPath filesep 'summary_' metadata{i,1} '.mat']);
    for j=1:length(sumResults)
        currRxn=model.rxns(sumResults(j,1));
        currRxn=strrep(currRxn,'extended_','');
        currRxn=strrep(currRxn,'[u]tr','');
        if contains(currRxn{1},'_IEX_')
            currRxnRow=find(strcmp(currRxn,ContributionFluxes(:,1)));
            if ~isempty(currRxnRow)
                ContributionFluxes{currRxnRow,1}=currRxn{1};
                ContributionFluxes{currRxnRow,i+8}=sumResults(j,2);
            else
                ContributionFluxes{cnt,1}=currRxn{1};
                % add the taxonomical information and metabolite produced
                [currStrainProd]=strsplit(currRxn{1},'_IEX_');
                ContributionFluxes{cnt,2}=currStrainProd{1,1};
                ContributionFluxes{cnt,3}=strrep(currStrainProd{1,2},'[u]tr','');
                for t=1:length(taxa)
                    taxCol=find(strcmp(taxa{t},infoFile(1,:)));
                    ContributionFluxes{cnt,t+3}=infoFile(find(strcmp(currStrainProd{1,1},infoFile(:,1))),taxCol);
                end
                ContributionFluxes{cnt,i+8}=sumResults(j,2);
                cnt=cnt+1;
            end
        end
    end
end
% replace empty entries with zeros to prevent formatting issues
% also remove positive values == uptake
for i=2:size(ContributionFluxes,1)
    for j=10:size(ContributionFluxes,2)
        if isempty(ContributionFluxes{i,j})
            ContributionFluxes{i,j}=0;
        end
        if ContributionFluxes{i,j}>0
            ContributionFluxes{i,j}=0;
        end
    end
end

% delete rows with very small/all zero values to reduce table size
cnt=1;
delArray=[];
for i=2:size(ContributionFluxes,1)
    sumRows=sum(abs(str2double(string(ContributionFluxes(i,10:end)))));
    if sumRows<0.000001
        delArray(cnt,1)=i;
        cnt=cnt+1;
    end
end
ContributionFluxes(delArray,:)=[];

ContributionFluxes(2:end,10:end)=num2cell(abs(cell2mat(ContributionFluxes(2:end,10:end))));
Taxonomy=ContributionFluxes(:,1:9);
Taxonomy{1,1}='Feature';
for i=1:size(Taxonomy,1)
    for j=1:size(Taxonomy,2)
        Taxonomy{i,j}=cellstr(Taxonomy{i,j});
        Taxonomy{i,j}=Taxonomy{i,j}{1,1};
    end
end
ContributionFluxes(:,2:9)=[];
cell2csv([contrPath filesep 'MicrobeContributions_Fluxes.csv'],ContributionFluxes);
cell2csv([contrPath filesep 'Taxonomy_MicrobeContributions_Fluxes.csv'],Taxonomy);
