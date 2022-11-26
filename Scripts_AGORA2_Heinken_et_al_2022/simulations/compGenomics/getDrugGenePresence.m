% Computes the presence of each drug gene in the ~5,000 analyzed genes.
% Generates the data underlying Figures 4b-d.

% get AGORA2 IDs
infoFile = readInputTableForPipeline('AGORA2_infoFile.xlsx');

spreadsheet=readInputTableForPipeline([spreadsheetFolder filesep 'M_Drugs.tsv']);

% get strains with comparative genomic data
genomeAnnotation = readtable('gapfilledGenomeAnnotation.txt', 'ReadVariableNames', false, 'Delimiter', 'tab','TreatAsEmpty',['UND. -60001','UND. -2011','UND. -62011']);
genomeAnnotation = table2cell(genomeAnnotation);
strains=unique(genomeAnnotation(:,1));

genePresence = {};
genePresence(2:length(strains)+1,1)=strains;

% delete unreconstructed strains
[C,IA]=setdiff(spreadsheet(:,1),infoFile(:,2),'stable');
spreadsheet(IA(2:end),:)=[];
cnt=2;

for i=2:size(spreadsheet,2)
    genePresence{1,cnt}=spreadsheet{1,i};
    genePresence(2:end,cnt)={'0'};
    for k=2:size(spreadsheet,1)
        findOrg=find(strcmp(infoFile(:,2),spreadsheet{k,1}));
        findInS=find(strcmp(genePresence(:,1),infoFile{findOrg,1}));
        if ~isempty(findInS)
            if ~isempty(spreadsheet{k,i})
                genePresence{findInS,cnt}='1';
            else
                genePresence{findInS,cnt}='0';
            end
        end
    end
    cnt=cnt+1;
end
cell2csv([propertiesFolder filesep 'AGORA2_DrugGenePresence.csv'],genePresence);
