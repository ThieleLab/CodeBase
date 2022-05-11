% Computes the presence of each gene in the ~5,000 analyzed genes

spreadsheetFolder=fileparts(which('B_Rib.tsv'));

% get AGORA2 IDs
infoFile = readtable('AGORA2_infoFile.xlsx', 'ReadVariableNames', false);
infoFile=table2cell(infoFile);

spreadsheet=readtable([spreadsheetFolder filesep 'M_Drugs.tsv'], 'ReadVariableNames', false,'FileType','text','delimiter','tab');
spreadsheet = table2cell(spreadsheet);
% delete unreconstructed strains
[C,IA]=setdiff(spreadsheet(:,1),infoFile(:,2),'stable');
spreadsheet(IA(2:end),:)=[];
for i=2:size(spreadsheet,1)
    spreadsheet{i,1}=infoFile{find(strcmp(infoFile(:,2),spreadsheet{i,1})),1};
end

taxLevels={'Phylum','Class','Order','Family','Genus','Species'};
for t=1:length(taxLevels)
    genePresence={taxLevels{t},'Total number'};
    getTax={};
    taxCol=find(strcmp(infoFile(1,:),taxLevels{t}));
    for k=2:size(spreadsheet,1)
        getTax{k,1}=infoFile{find(strcmp(infoFile(:,1),spreadsheet{k,1})),taxCol};
    end
    [taxa, ~, J]=unique(getTax(2:end,1));
    taxnum = histc(J, 1:numel(taxa));
    for j=1:length(taxa)
        genePresence{j+1,1}=taxa{j};
        genePresence{j+1,2}=taxnum(j);
    end
    for i=2:size(spreadsheet,2)
        genePresence{1,i+1}=spreadsheet{1,i};
        genePresence(2:end,i+1)={0};
        % get all strains
        genes=spreadsheet(:,1);
        genes(:,2)=spreadsheet(:,i);
        cnt=1;
        delArray=[];
        for k=2:size(genes,1)
            if isempty(genes{k,2})
                delArray(cnt,1)=k;
                cnt=cnt+1;
            end
        end
        genes(delArray,:)=[];
        % get taxonomic distribution
        taxCol=find(strcmp(infoFile(1,:),taxLevels{t}));
        for k=2:size(genes,1)
            genes{k,3}=infoFile{find(strcmp(infoFile(:,1),genes{k,1})),taxCol};
        end
        [taxa, ~, J]=unique(genes(2:end,3));
        taxnum = histc(J, 1:numel(taxa));
        for j=1:length(taxa)
            findRow=find(strcmp(genePresence(:,1),taxa{j}));
            genePresence{findRow,i+1}=taxnum(j);
        end
    end
    genePresence(find(strncmp(genePresence(:,1),'unclassified',length('unclassified'))),:)=[];
    genePresencePercentage=genePresence;
    genePresence=cell2table(genePresence);
    writetable(genePresence,[resultsFolder 'DrugGenePresence_' taxLevels{t}],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
    
    for i=3:size(genePresencePercentage,2)
        for j=2:size(genePresencePercentage,1)
            genePresencePercentage{j,i}=genePresencePercentage{j,i}/genePresencePercentage{j,2};
        end
    end
    genePresencePercentage(:,2)=[];
    genePresencePercentage=cell2table(genePresencePercentage);
    writetable(genePresencePercentage,[propertiesFolder 'DrugGenePresence_Percentage_' taxLevels{t}],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
end
