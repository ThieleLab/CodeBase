% exports all computed data points in comparison for Table S5

cd([rootDir filesep 'Comparison_other_GEMs'])

resources={
    'NJC19','Uptake_Lim','Secretion_Lim'
    'BacDive_Metabolites','Uptake_BacDive','Secretion_BacDive'
    'BacDive_Enzymes','Enzymes_BacDive',''
    'Madin','Uptake_Madin',''
    };

for i=1:size(resources,1)
    load([rootDir filesep 'Comparison_other_GEMs' filesep 'AGORA2' filesep 'Comparison_Findings.mat']);
    
    Table={'Organism','Exchange','Direction','AGORA2','KBase','BiGG','CarveMe','gapseq','MAGMA'};
    
    orgs={};
    for j=2:size(resources,2)
        if ~strcmp(resources{i,j},'')
        orgs=union(orgs,unique(findings.(resources{i,j}).('AGORA2')(:,1)));
        end
    end
    
    cnt=2;
    
    for k=1:length(orgs)
        for m=2:size(resources,2)
            if ~strcmp(resources{i,m},'')
                findCases=find(strcmp(findings.(resources{i,m}).('AGORA2')(:,1),orgs{k}));
                for l=1:length(findCases)
                    Table{cnt,1}=findings.(resources{i,m}).('AGORA2'){findCases(l),1};
                    Table{cnt,2}=findings.(resources{i,m}).('AGORA2'){findCases(l),3};
                    if m==2
                        Table{cnt,3}='Uptake';
                    elseif m==3
                        Table{cnt,3}='Secretion';
                    end
                    Table{cnt,4}=findings.(resources{i,m}).('AGORA2'){findCases(l),2};
                    Table{cnt,5}=findings.(resources{i,m}).('KBase'){findCases(l),2};
                    cnt=cnt+1;
                end
            end
        end
    end
    
    % find for the other resources if possible
    for j=6:size(Table,2)
        load([rootDir filesep 'Comparison_other_GEMs' filesep Table{1,j} filesep 'Comparison_Findings.mat']);
        if size(findings.(resources{i,2}).(Table{1,j}),1)>0
            for k=2:size(Table,1)
                if strcmp(Table{k,3},'Uptake')
                    findOrg=find(strcmp(findings.(resources{i,2}).(Table{1,j})(:,1),Table{k,1}));
                    findMet=find(strcmp(findings.(resources{i,2}).(Table{1,j})(:,3),Table{k,2}));
                    if ~isempty(intersect(findOrg,findMet))
                        Table{k,j}=findings.(resources{i,2}).(Table{1,j}){intersect(findOrg,findMet),2};
                    else
                        Table{k,j}='N/A';
                    end
                elseif strcmp(Table{k,3},'Secretion')
                    findOrg=find(strcmp(findings.(resources{i,3}).(Table{1,j})(:,1),Table{k,1}));
                    findMet=find(strcmp(findings.(resources{i,3}).(Table{1,j})(:,3),Table{k,2}));
                    if ~isempty(intersect(findOrg,findMet))
                        Table{k,j}=findings.(resources{i,3}).(Table{1,j}){intersect(findOrg,findMet),2};
                    else
                        Table{k,j}='N/A';
                    end
                end
            end
        end
    end
    cell2csv([pwd filesep 'Results' filesep resources{i,1} '_summary.csv'],Table);
end
