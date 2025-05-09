function bestCorrPerSizeList = addCorrForSummedFluxLinkedTaxa(corrTable,bestCorrPerSizeList, rxnsOfInterest)

% Remove sum of taxa row
fluxSetCorr = corrTable(~matches(corrTable.Row, 'Sum of taxa'),:);

% Add correlations of the full set of flux-associated taxa
fluxSetCorr = fluxSetCorr(matches(fluxSetCorr.Row, 'Flux-associated taxa'),:);

% Get the number of flux associated microbes
numMicrobes = sum(corrTable{:,:}~=0,1)-1;

for i=1:length(bestCorrPerSizeList)
    bestCorrPerSize = bestCorrPerSizeList{i};

    % Get the first row (correlation results for the single best
    % correlating microbe)
    firstRow = bestCorrPerSize(1,:);
    firstRow.Microbe_1 = ""; % Remove associated microbe
    % Add cluster size and RHO for correlation with all combined relative
    % abunbdances o flux-associated microbes
    firstRow{:,{'Cluster_size','Spearman_RHO'}} = [numMicrobes(i) fluxSetCorr.(rxnsOfInterest{i})]; 

    % Extend bestCorrPerSize
    bestCorrPerSizeList{i} = [bestCorrPerSize; firstRow];

    % Add the total correlation coefficients
    % warning('off')
    % bestCorrPerSize.(rxnsOfInterest{i})(end+1) = fluxSetCorr.(rxnsOfInterest{i});
    % % Set rho_improvement of all taxa to nan
    % bestCorrPerSize.rho_improvement(end) = nan;
    % warning('on')
    % 
    % % Label new row
    % bestCorrPerSize(end,1) = {'Flux-associated taxa'};
    % 
    % % Add the total number of microbes
    % bestCorrPerSize.Size(end) = numMicrobes(i);
    % 
    % % Update table
    % bestCorrPerSizeList{i} = bestCorrPerSize;
end
end