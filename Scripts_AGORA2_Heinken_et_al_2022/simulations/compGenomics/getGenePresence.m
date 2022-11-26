% Computes the presence of each gene in the >5,000 analyzed strains
% Generates Table S3c

dInfo = dir(spreadsheetFolder);
fileList={dInfo.name};
fileList=fileList';
fileList(find(strcmp(fileList(:,1),'.')),:)=[];
fileList(find(strcmp(fileList(:,1),'..')),:)=[];

% get AGORA2 IDs
infoFile = readInputTableForPipeline('AGORA2_infoFile.xlsx');

genePresence={};
geneSubsystem={'Gene','Subsystem'};

% get strains with comparative genomic data
genomeAnnotation = readtable('gapfilledGenomeAnnotation.txt', 'ReadVariableNames', false, 'Delimiter', 'tab','TreatAsEmpty',['UND. -60001','UND. -2011','UND. -62011']);
genomeAnnotation = table2cell(genomeAnnotation);
strains=unique(genomeAnnotation(:,1));
genePresence(2:length(strains)+1,1)=strains;

cnt=2;
for i=1:length(fileList)
    i
    spreadsheet=readInputTableForPipeline([spreadsheetFolder filesep fileList{i}]);
    % delete unreconstructed strains
    [C,IA]=setdiff(spreadsheet(:,1),infoFile(:,2),'stable');
    spreadsheet(IA(2:end),:)=[];
    for j=2:size(spreadsheet,2)
        genePresence{1,cnt}=spreadsheet{1,j};
        genePresence(2:end,cnt)={'0'};
        for k=2:size(spreadsheet,1)
            findOrg=find(strcmp(infoFile(:,2),spreadsheet{k,1}));
            findInS=find(strcmp(genePresence(:,1),infoFile{findOrg,1}));
            if ~isempty(findInS)
                if ~isempty(spreadsheet{k,j})
                    genePresence{findInS,cnt}='1';
                else
                    genePresence{findInS,cnt}='0';
                end
            end
        end
        cnt=cnt+1;
    end
end

% find and remove duplicate columns
genePresence{1,1}='Strain';
[C,IA,IB]  = unique(genePresence(1,:));
repeatedStr = C(histc(IB,1:max(IB))>1);
duplicateStrains=repeatedStr;
delArray=[];
cnt=1;
if ~isempty(repeatedStr)
    for i=1:length(repeatedStr)
        countdata=[];
        findintable=find(strcmp(genePresence(1,:),repeatedStr{i}));
        % delete the one with less or no data
        for j=1:length(findintable)
            countdata(j)=sum(str2double(genePresence(2:end,findintable(j))));
        end
        % entry with most data points will not be deleted (or first one is
        % kept if all are zero)
        [M,I]=max(countdata);
        findintable(I)=[];
       for j=1:length(findintable)
           delArray(cnt,1)=findintable(j);
           cnt=cnt+1;
       end
    end
end
genePresence(:,delArray)=[];

cell2csv([propertiesFolder filesep 'AGORA2_GenePresence.csv'],genePresence);

