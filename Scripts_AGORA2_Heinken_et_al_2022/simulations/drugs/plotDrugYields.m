
resultsFolder = [rootDir filesep 'ComputedDrugFluxes' filesep];

% plot drug yields

yields={
    'ATP' 'atp'
    'CO2' 'co2'
    'Pyruvate' 'pyr'
    'NH4' 'nh4'
    };

cols = zeros(size(yields,2),3);

f=figure;
for i=1:length(yields)
    yield = readInputTableForPipeline([resultsFolder 'AGORA2_drugYields_' yields{i,2} '.txt']);
    yield(1,:)=strrep(yield(1,:),'5-Aminosalicylic Acid (Mesalamine)','5-aminosalicylic acid');
    cnt=1;
    data=[];
    labels={};
    for j=2:size(yield,2)
        for k=2:size(yield,1)
            if j==2 || yield{k,j} > 0.000001
                data(cnt,1)=yield{k,j};
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
    set(gca, 'FontSize', 12)
%     box on
    xtickangle(45)
    set(gca,'TickLabelInterpreter','none')
    title([yields{i,1} ' yield' ', n = ' num2str(size(yield,1)-1)])
end
f.Renderer='painters';
