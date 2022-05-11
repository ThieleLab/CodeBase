                                                                                                                                                            
% Get averages and standard deviations for basic reconstruction features
% Draft and curated reconstructions are compared

% Versions
versions = {'AGORA_2_0_Curated','AGORA_2_0_Draft'};
files = {[propertiesFolder filesep 'All_statistics_Curated.csv'],[propertiesFolder filesep 'All_statistics_Draft.csv']};

for i=1:2
    Averages{1,i+1} = versions{i};
    stats = readtable(files{i}, 'ReadVariableNames', false);
    stats = table2cell(stats);
    for j=2:size(stats,2)
        Averages{j,1} = stats{1,j};
        if strncmp(stats{1,j},'Growth',6)
            Averages{j,i+1} = num2str(sum(str2double(stats(2:end,j))> 0.000001));
        else
            av = mean(str2double(stats(2:end,j)));
            s = std(str2double(stats(2:end,j)));
            Averages{j,i+1} = [num2str(round(av,2)) ' +/- ' num2str(round(s,2))];
        end
    end
end

