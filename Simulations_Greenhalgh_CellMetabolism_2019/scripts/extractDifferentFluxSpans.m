%% new tables 25.04.2018

%% all host reactions that have >20% difference in flux span in at least one condition
load('FluxSpans_DMEM_LGG.mat');
FluxSpans_DMEM=FluxSpans;
load('FluxSpans_SIEM_LGG.mat');
FluxSpans_SIEM=FluxSpans;
Reactions={};
Reactions{1,1}='VMH_ID';
Reactions{1,2}='Description';
Reactions{1,3}='Subsystem';
Reactions{1,4}='Caco2_FF';
Reactions{1,5}='Caco2_LGG_FF';
Reactions{1,6}='Caco2_HF';
Reactions{1,7}='Caco2_LGG_HF';
Reactions{1,8}='Model_RxnID';
cnt=2;
for i=1:length(FluxSpans_DMEM)
    if strncmp(FluxSpans_DMEM{i,1},'Host_',5)
        % get all values
        vals(1,1)=FluxSpans_DMEM{i,9};
        vals(1,2)=FluxSpans_DMEM{i,10};
        vals(1,3)=FluxSpans_SIEM{i,9};
        vals(1,4)=FluxSpans_SIEM{i,10};
        % only those with somehow noticable flux span
        diff=false;
        if vals(1,:) >0.01
            % find the ones with at least 20% difference between any two conditions
            for j=1:4
                for k=2:4
                    if vals(k)>0 && vals(j)/vals(k)>1.2 || vals(j)>0 && vals(k)/vals(j)>1.2
                        diff=true;
                    end
                end
            end
        end
        if diff==true
            Reactions{cnt,1}=strrep(FluxSpans_DMEM{i,1},'Host_','');
            Reactions{cnt,2}=FluxSpans_DMEM{i,3};
            Reactions{cnt,3}=FluxSpans_DMEM{i,5};
            Reactions{cnt,4}=FluxSpans_DMEM{i,9};
            Reactions{cnt,5}=FluxSpans_DMEM{i,10};
            Reactions{cnt,6}=FluxSpans_SIEM{i,9};
            Reactions{cnt,7}=FluxSpans_SIEM{i,10};
            Reactions{cnt,8}=FluxSpans_DMEM{i,1};
            cnt=cnt+1;
        end
    end
end
save('FluxSpans_all_Different_Reactions_Host','Reactions');

%%
load('MaxFluxes_DMEM_LGG.mat');
MaxFluxes_DMEM=MaxFluxes;
load('MaxFluxes_SIEM_LGG.mat');
MaxFluxes_SIEM=MaxFluxes;
Reactions={};
Reactions{1,1}='VMH_ID';
Reactions{1,2}='Description';
Reactions{1,3}='Subsystem';
Reactions{1,4}='Caco2_FF';
Reactions{1,5}='Caco2_LGG_FF';
Reactions{1,6}='Caco2_HF';
Reactions{1,7}='Caco2_LGG_HF';
Reactions{1,8}='Model_RxnID';
cnt=2;
for i=1:length(MaxFluxes_DMEM)
    if strncmp(MaxFluxes_DMEM{i,1},'Host_',5)
        % get all values
        vals(1,1)=MaxFluxes_DMEM{i,9};
        vals(1,2)=MaxFluxes_DMEM{i,10};
        vals(1,3)=MaxFluxes_SIEM{i,9};
        vals(1,4)=MaxFluxes_SIEM{i,10};
        % only those with somehow noticable flux span
        diff=false;
        if vals(1,:) >0.01
            % find the ones with at least 20% difference between any two conditions
            for j=1:4
                for k=2:4
                    if vals(k)>0 && vals(j)/vals(k)>1.2 || vals(j)>0 && vals(k)/vals(j)>1.2
                        diff=true;
                    end
                end
            end
        end
        if diff==true
            Reactions{cnt,1}=strrep(MaxFluxes_DMEM{i,1},'Host_','');
            Reactions{cnt,2}=MaxFluxes_DMEM{i,3};
            Reactions{cnt,3}=MaxFluxes_DMEM{i,5};
            Reactions{cnt,4}=MaxFluxes_DMEM{i,9};
            Reactions{cnt,5}=MaxFluxes_DMEM{i,10};
            Reactions{cnt,6}=MaxFluxes_SIEM{i,9};
            Reactions{cnt,7}=MaxFluxes_SIEM{i,10};
            Reactions{cnt,8}=MaxFluxes_DMEM{i,1};
            cnt=cnt+1;
        end
    end
end
save('MaxFluxes_all_Different_Reactions_Host','Reactions');

%% all LGG reactions that have >20% difference in flux span in at least one condition
load('FluxSpans_DMEM_LGG.mat');
FluxSpans_DMEM=FluxSpans;
load('FluxSpans_SIEM_LGG.mat');
FluxSpans_SIEM=FluxSpans;
Reactions={};
Reactions{1,1}='VMH_ID';
Reactions{1,2}='Description';
Reactions{1,3}='Subsystem';
Reactions{1,4}='LGG_FF';
Reactions{1,5}='Caco2_LGG_FF';
Reactions{1,6}='LGG_HF';
Reactions{1,7}='Caco2_LGG_HF';
Reactions{1,8}='Model_RxnID';
cnt=2;
for i=1:length(FluxSpans_DMEM)
    if strncmp(FluxSpans_DMEM{i,1},'Lactobacillus_rhamnosus_GG_ATCC_53103_',length('Lactobacillus_rhamnosus_GG_ATCC_53103_'))
        % get all values
        vals(1,1)=FluxSpans_DMEM{i,7};
        vals(1,2)=FluxSpans_DMEM{i,10};
        vals(1,3)=FluxSpans_SIEM{i,7};
        vals(1,4)=FluxSpans_SIEM{i,10};
        % only those with somehow noticable flux span
        diff=false;
        if vals(1,:) >0.01
            % find the ones with at least 20% difference between any two conditions
            for j=1:4
                for k=2:4
                    if vals(k)>0 && vals(j)/vals(k)>1.2 || vals(j)>0 && vals(k)/vals(j)>1.2
                        diff=true;
                    end
                end
            end
        end
        if diff==true
            Reactions{cnt,1}=strrep(FluxSpans_DMEM{i,1},'Lactobacillus_rhamnosus_GG_ATCC_53103_','');
            Reactions{cnt,2}=FluxSpans_DMEM{i,3};
            Reactions{cnt,3}=FluxSpans_DMEM{i,5};
            Reactions{cnt,4}=FluxSpans_DMEM{i,7};
            Reactions{cnt,5}=FluxSpans_DMEM{i,10};
            Reactions{cnt,6}=FluxSpans_SIEM{i,7};
            Reactions{cnt,7}=FluxSpans_SIEM{i,10};
            Reactions{cnt,8}=FluxSpans_DMEM{i,1};
            cnt=cnt+1;
        end
    end
end
save('FluxSpans_all_Different_Reactions_LGG','Reactions');

%%
load('MaxFluxes_DMEM_LGG.mat');
MaxFluxes_DMEM=MaxFluxes;
load('MaxFluxes_SIEM_LGG.mat');
MaxFluxes_SIEM=MaxFluxes;
Reactions={};
Reactions{1,1}='VMH_ID';
Reactions{1,2}='Description';
Reactions{1,3}='Subsystem';
Reactions{1,4}='LGG_FF';
Reactions{1,5}='LGG_Caco2_FF';
Reactions{1,6}='LGG_HF';
Reactions{1,7}='LGG_Caco_2HF';
Reactions{1,8}='Model_RxnID';
cnt=2;
for i=1:length(MaxFluxes_DMEM)
    if strncmp(MaxFluxes_DMEM{i,1},'Lactobacillus_rhamnosus_GG_ATCC_53103_',length('Lactobacillus_rhamnosus_GG_ATCC_53103_'))
        % get all values
        vals(1,1)=MaxFluxes_DMEM{i,7};
        vals(1,2)=MaxFluxes_DMEM{i,10};
        vals(1,3)=MaxFluxes_SIEM{i,7};
        vals(1,4)=MaxFluxes_SIEM{i,10};
        % only those with somehow noticable flux span
        diff=false;
        if vals(1,:) >0.01
            % find the ones with at least 20% difference between any two conditions
            for j=1:4
                for k=2:4
                    if vals(k)>0 && vals(j)/vals(k)>1.2 || vals(j)>0 && vals(k)/vals(j)>1.2
                        diff=true;
                    end
                end
            end
        end
        if diff==true
            Reactions{cnt,1}=strrep(MaxFluxes_DMEM{i,1},'Lactobacillus_rhamnosus_GG_ATCC_53103_','');
            Reactions{cnt,2}=MaxFluxes_DMEM{i,3};
            Reactions{cnt,3}=MaxFluxes_DMEM{i,5};
            Reactions{cnt,4}=MaxFluxes_DMEM{i,7};
            Reactions{cnt,5}=MaxFluxes_DMEM{i,10};
            Reactions{cnt,6}=MaxFluxes_SIEM{i,7};
            Reactions{cnt,7}=MaxFluxes_SIEM{i,10};
            Reactions{cnt,8}=MaxFluxes_DMEM{i,1};
            cnt=cnt+1;
        end
    end
end
save('MaxFluxes_all_Different_Reactions_LGG','Reactions');

%% count the statistics
load('FluxSpans_DMEM_LGG.mat');
FluxSpans_DMEM=FluxSpans;
load('FluxSpans_SIEM_LGG.mat');
FluxSpans_SIEM=FluxSpans;

subs=unique(cellstr(string(FluxSpans_DMEM(2:end,5))));
Table={};
Table{1,2}='Caco2_LGG_FF_vs_Caco2_FF';
Table{1,3}='Caco2_FF_vs_Caco2_LGG_FF';
Table{1,4}='Caco2_LGG_HF_vs_Caco2_HF';
Table{1,5}='Caco2_HF_vs_Caco2_LGG_HF';
Table{1,6}='Caco2_LGG_FF_vs_Caco2_LGG_HF';
Table{1,7}='Caco2_LGG_HF_vs_Caco2_LGG_FF';
Table{1,8}='Caco2_FF_vs_Caco2_HF';
Table{1,9}='Caco2_HF_vs_Caco2_FF';
Table{1,10}='LGG_Caco2_FF_vs_LGG_FF';
Table{1,11}='LGG_FF_vs_LGG_Caco2_FF';
Table{1,12}='LGG_Caco2_HF_vs_LGG_HF';
Table{1,13}='LGG_HF_vs_LGG_Caco2_HF';
Table{1,14}='LGG_Caco2_FF_vs_LGG_Caco2_HF';
Table{1,15}='LGG_Caco2_HF_vs_LGG_Caco2_FF';
Table{1,16}='LGG_FF_vs_LGG_HF';
Table{1,17}='LGG_HF_vs_LGG_FF';
for i=1:length(subs)
    Table{i+1,1}=subs{i};
    for j=2:size(Table,2)
        Table{i+1,j}=0;
    end
end
for i=1:length(FluxSpans_DMEM)
    % first host
    if strncmp(FluxSpans_DMEM{i,1},'Host_',5)
        % get all values
        val_Single_DMEM=FluxSpans_DMEM{i,9};
        val_Co_DMEM=FluxSpans_DMEM{i,10};
        val_Single_SIEM=FluxSpans_SIEM{i,9};
        val_Co_SIEM=FluxSpans_SIEM{i,10};
        % find the ones with at least 20% difference between any two conditions
        if val_Single_DMEM >0.01 || val_Co_DMEM >0.01 || val_Single_SIEM >0.01 || val_Co_SIEM >0.01
            % find the subsystem
            sub=find(strcmp(Table(:,1),FluxSpans_DMEM{i,5}));
            % Comparisons
            if  val_Co_DMEM/val_Single_DMEM >1.2 || val_Co_DMEM>0 && val_Single_DMEM<0.00000001
                Table{sub,2}=Table{sub,2}+1;
            end
            if  val_Single_DMEM/val_Co_DMEM >1.2 || val_Single_DMEM>0 && val_Co_DMEM<0.00000001
                Table{sub,3}=Table{sub,3}+1;
            end
            if  val_Co_SIEM/val_Single_SIEM >1.2 || val_Co_SIEM>0 && val_Single_SIEM<0.00000001
                Table{sub,4}=Table{sub,4}+1;
            end
            if  val_Single_SIEM/val_Co_SIEM >1.2 || val_Single_SIEM>0 && val_Co_SIEM<0.00000001
                Table{sub,5}=Table{sub,5}+1;
            end
            if  val_Co_DMEM/val_Co_SIEM >1.2 || val_Co_DMEM>0 && val_Co_SIEM<0.00000001
                Table{sub,6}=Table{sub,6}+1;
            end
            if  val_Co_SIEM/val_Co_DMEM >1.2 || val_Co_SIEM>0 && val_Co_DMEM<0.00000001
                Table{sub,7}=Table{sub,7}+1;
            end
            if  val_Single_DMEM/val_Single_SIEM >1.2 || val_Single_DMEM>0 && val_Single_SIEM<0.00000001
                Table{sub,8}=Table{sub,8}+1;
            end
            if  val_Single_SIEM/val_Single_DMEM >1.2 || val_Single_SIEM>0 && val_Single_DMEM<0.00000001
                Table{sub,9}=Table{sub,9}+1;
            end
        end
    end
    % then LGG
    if strncmp(FluxSpans_DMEM{i,1},'Lactobacillus_rhamnosus_GG_ATCC_53103_',length('Lactobacillus_rhamnosus_GG_ATCC_53103_'))
        % get all values
        val_Single_DMEM=FluxSpans_DMEM{i,7};
        val_Co_DMEM=FluxSpans_DMEM{i,10};
        val_Single_SIEM=FluxSpans_SIEM{i,7};
        val_Co_SIEM=FluxSpans_SIEM{i,10};
        % find the ones with at least 20% difference between any two conditions
        if val_Single_DMEM >0.01 || val_Co_DMEM >0.01 || val_Single_SIEM >0.01 || val_Co_SIEM >0.01
            % find the subsystem
            sub=find(strcmp(Table(:,1),FluxSpans_DMEM{i,5}));
            % Comparisons
            if  val_Co_DMEM/val_Single_DMEM >1.2 || val_Co_DMEM>0 && val_Single_DMEM<0.00000001
                Table{sub,10}=Table{sub,10}+1;
            end
            if  val_Single_DMEM/val_Co_DMEM >1.2 || val_Single_DMEM>0 && val_Co_DMEM<0.00000001
                Table{sub,11}=Table{sub,11}+1;
            end
            if  val_Co_SIEM/val_Single_SIEM >1.2 || val_Co_SIEM>0 && val_Single_SIEM<0.00000001
                Table{sub,12}=Table{sub,12}+1;
            end
            if  val_Single_SIEM/val_Co_SIEM >1.2 || val_Single_SIEM>0 && val_Co_SIEM<0.00000001
                Table{sub,13}=Table{sub,13}+1;
            end
            if  val_Co_DMEM/val_Co_SIEM >1.2 || val_Co_DMEM>0 && val_Co_SIEM<0.00000001
                Table{sub,14}=Table{sub,14}+1;
            end
            if  val_Co_SIEM/val_Co_DMEM >1.2 || val_Co_SIEM>0 && val_Co_DMEM<0.00000001
                Table{sub,15}=Table{sub,15}+1;
            end
            if  val_Single_DMEM/val_Single_SIEM >1.2 || val_Single_DMEM>0 && val_Single_SIEM<0.00000001
                Table{sub,16}=Table{sub,16}+1;
            end
            if  val_Single_SIEM/val_Single_DMEM >1.2 || val_Single_SIEM>0 && val_Single_DMEM<0.00000001
                Table{sub,17}=Table{sub,17}+1;
            end
        end
    end
end

% delete zero rows/columns
cnt=1;
for i=2:size(Table,1)
    summed=sum(cell2mat(Table(i,2:end)));
    if summed==0
        delRow(cnt,1)=i;
        cnt=cnt+1;
    end
end
Table(delRow,:)=[];

cnt=1;
for i=2:size(Table,2)
    summed=sum(cell2mat(Table(2:end,i)));
    if summed==0
        delCol(1,cnt)=i;
        cnt=cnt+1;
    end
end
Table(:,delCol)=[];

