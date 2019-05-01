%% Extract reactions with differential flux spans in LGG and Caco-2 (host). Related to Figures 2A, 4C, 5B-C

%% all host reactions that have >20% difference in flux span in at least one condition
load('FluxSpans_REF_LGG.mat');
FluxSpans_REF=FluxSpans;
load('FluxSpans_HF_LGG.mat');
FluxSpans_HF=FluxSpans;
Reactions={};
Reactions{1,1}='Model_RxnID';
Reactions{1,2}='VMH_ID';
Reactions{1,3}='Description';
Reactions{1,4}='Subsystem';
Reactions{1,5}='Caco2_FF';
Reactions{1,6}='Caco2_LGG_FF';
Reactions{1,7}='Caco2_HF';
Reactions{1,8}='Caco2_LGG_HF';
cnt=2;
for i=1:length(FluxSpans_REF)
    if strncmp(FluxSpans_REF{i,1},'Host_',5)
        % get all values
        vals(1,1)=FluxSpans_REF{i,8};
        vals(1,2)=FluxSpans_REF{i,9};
        vals(1,3)=FluxSpans_HF{i,8};
        vals(1,4)=FluxSpans_HF{i,9};
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
            Reactions{cnt,1}=FluxSpans_REF{i,1};
            Reactions{cnt,2}=strrep(FluxSpans_REF{i,1},'Host_','');
            Reactions{cnt,3}=FluxSpans_REF{i,3};
            Reactions{cnt,4}=FluxSpans_REF{i,5};
            Reactions{cnt,5}=FluxSpans_REF{i,8};
            Reactions{cnt,6}=FluxSpans_REF{i,9};
            Reactions{cnt,7}=FluxSpans_HF{i,8};
            Reactions{cnt,8}=FluxSpans_HF{i,9};
            cnt=cnt+1;
        end
    end
end
save('FluxSpans_all_Different_Reactions_Host','Reactions');

%% all LGG reactions that have >20% difference in flux span in at least one condition
load('FluxSpans_REF_LGG.mat');
FluxSpans_REF=FluxSpans;
load('FluxSpans_HF_LGG.mat');
FluxSpans_HF=FluxSpans;
Reactions={};
Reactions{1,1}='Model_RxnID';
Reactions{1,2}='VMH_ID';
Reactions{1,3}='Description';
Reactions{1,4}='Subsystem';
Reactions{1,5}='LGG_FF';
Reactions{1,6}='Caco2_LGG_FF';
Reactions{1,7}='LGG_HF';
Reactions{1,8}='Caco2_LGG_HF';
cnt=2;
for i=1:length(FluxSpans_REF)
    if strncmp(FluxSpans_REF{i,1},'Lactobacillus_rhamnosus_GG_ATCC_53103_',length('Lactobacillus_rhamnosus_GG_ATCC_53103_'))
        % get all values
        vals(1,1)=FluxSpans_REF{i,7};
        vals(1,2)=FluxSpans_REF{i,9};
        vals(1,3)=FluxSpans_HF{i,7};
        vals(1,4)=FluxSpans_HF{i,9};
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
            Reactions{cnt,1}=FluxSpans_REF{i,1};
            Reactions{cnt,2}=strrep(FluxSpans_REF{i,1},'Lactobacillus_rhamnosus_GG_ATCC_53103_','');
            Reactions{cnt,3}=FluxSpans_REF{i,3};
            Reactions{cnt,4}=FluxSpans_REF{i,5};
            Reactions{cnt,5}=FluxSpans_REF{i,7};
            Reactions{cnt,6}=FluxSpans_REF{i,9};
            Reactions{cnt,7}=FluxSpans_HF{i,7};
            Reactions{cnt,8}=FluxSpans_HF{i,9};
            cnt=cnt+1;
        end
    end
end
save('FluxSpans_all_Different_Reactions_LGG','Reactions');

