
resultsFolder = [rootDir filesep 'ComputeDrugReactions' filesep];
plotsFolder = [rootDir filesep 'ComputeDrugReactions' filesep 'Figures' filesep];

% plot drug yields

yields={
    'ATP' 'atp'
    'CO2' 'co2'
    'Pyruvate' 'pyr'
    'NH4' 'nh4'
    };

% plot any drugs that have a yield above zero
for i=1:length(yields)
    yield = readtable([resultsFolder 'AGORA2_drugYields_' yields{i,2} '.txt'], 'ReadVariableNames', false);
    yield = table2cell(yield);
    yield(1,:)=strrep(yield(1,:),'5-aminosalicylic acid (Mesalamine)','5-aminosalicylic acid');
    cnt=1;
    data=[];
    labels={};
    cols=[];
    for j=2:size(yield,2)
        for k=2:size(yield,1)
            if j==2 || str2double(yield{k,j}) > 0.000001
                data(cnt,1)=str2double(yield{k,j});
                labels{cnt,1}=strrep(yield{1,j},['_' yields{i,2}],'');
                cnt=cnt+1;
            end
        end
        cols(j,:)=[rand rand rand];
    end
    f=figure;
    boxplot(data,labels,'boxstyle','filled','colors',cols)
    ylim([0 max(data)+1])
    ylabel('Yield from 1 mmol*g dry weight-1*hr-1 drug')
    set(gca, 'FontSize', 10)
    box on
    xtickangle(45)
    set(gca,'TickLabelInterpreter','none')
    title([yields{i,1} ' yield' ', n = ' num2str(size(yield,1)-1)])
    f.Renderer='painters';
    print([plotsFolder yields{i,1} '_yield'],'-dpng','-r300')
end

% all in one figure

yields={
    'ATP' 'atp'
    'CO2' 'co2'
    'Pyruvate' 'pyr'
    'NH4' 'nh4'
    };

f=figure;
for i=1:length(yields)
    yield = readtable([resultsFolder 'AGORA2_drugYields_' yields{i,2} '.txt'], 'ReadVariableNames', false);
    yield = table2cell(yield);
    yield(1,:)=strrep(yield(1,:),'5-Aminosalicylic Acid (Mesalamine)','5-aminosalicylic acid');
    cnt=1;
    data=[];
    labels={};
    for j=2:size(yield,2)
        for k=2:size(yield,1)
            if j==2 || str2double(yield{k,j}) > 0.000001
                data(cnt,1)=str2double(yield{k,j});
                labels{cnt,1}=strrep(yield{1,j},['_' yields{i,2}],'');
                cnt=cnt+1;
            end
        end
        cols(j,:)=[rand rand rand];
    end
    subplot(2,2,i)
    boxplot(data,labels,'boxstyle','filled','colors',cols)
    ylim([0 max(data)+1])
    ylabel('mmol*g dry weight-1*hr-1')
    set(gca, 'FontSize', 10)
    box on
    xtickangle(45)
    set(gca,'TickLabelInterpreter','none')
    title([yields{i,1} ' yield' ', n = ' num2str(size(yield,1)-1)])
end
f.Renderer='painters';
print([plotsFolder 'Drug_yields'],'-dpng','-r600')

