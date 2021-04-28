% Extract strain contribution profiles for each metabolite

Fluxes=table2cell(readtable([contrPath filesep 'MicrobeContributions_Fluxes.csv'], 'ReadVariableNames', false));
Taxonomy=table2cell(readtable([contrPath filesep 'Taxonomy_MicrobeContributions_Fluxes.csv'], 'ReadVariableNames', false));

%% separate files for each metabolite
sheetnames={};
sheetcnt=1;
for i=1:length(metaboliteInfo)
    FluxesTemp=Fluxes;
    TaxonomyTemp=Taxonomy;
    % first extract all rows with that metabolite
    cnt=1;
    delArray=[];
    for j=2:size(FluxesTemp,1)
        if ~strcmp(TaxonomyTemp(j,3),metaboliteInfo(i,2))
            delArray(cnt,1)=j;
            cnt=cnt+1;
        end
    end
    FluxesTemp(delArray,:)=[];
    TaxonomyTemp(delArray,:)=[];
    FluxesTemp(2:end,1)=strrep(FluxesTemp(2:end,1),strcat('_IEX_',metaboliteInfo{i,2}),'');
    % Only proceed if contributions not zero or too small, also remove very
    % small individual contributions
    sums=sum(str2double(FluxesTemp(2:end,2:end)'));
    if sum(sums) > 0.01
        % remove all rows with small values
        cnt=1;
        delArray=[];
        % find largest row
        for j=2:size(FluxesTemp,1)
            % remove all that contribute less than 1% of largest row summed up
            if sum(str2double(FluxesTemp(j,2:end))) <max(sums)*0.1
                delArray(cnt,1)=j;
                cnt=cnt+1;
            end
        end
        FluxesTemp(delArray,:)=[];
        TaxonomyTemp(delArray,:)=[];
        % Keep only metabolites with more than five strains
        % contributing-hard to plot otherwise
        if size(FluxesTemp,1) > 5
            filename=strcat('Secretion_',metaboliteInfo(i,2),'.csv');
            sheetnames{sheetcnt,1}=sheetcnt;
            sheetnames{sheetcnt,2}=filename;
            cell2csv([profilePath filesep filename{1}],FluxesTemp);
            filename=strcat('Taxonomy_secretion_',metaboliteInfo(i,2),'.csv');
            sheetnames{sheetcnt,3}=filename;
            sheetnames{sheetcnt,4}=metaboliteInfo(i,1);
            sheetnames{sheetcnt,5}=metaboliteInfo(i,2);
            cell2csv([profilePath filesep filename{1}],TaxonomyTemp);
            sheetcnt=sheetcnt+1;
        end
    end
end
for i=1:size(sheetnames,1)
    for j=2:size(sheetnames,2)
    sheetnames(i,j)=sheetnames{i,j};
    end
end
cell2csv([profilePath filesep 'sheetnames_secretion_metabolites.csv'],sheetnames);

clear Fluxes
clear Taxonomy
