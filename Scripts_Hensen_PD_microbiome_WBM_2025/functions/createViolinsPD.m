function fig = createViolinsPD(regressionResults,preparedInputTable,preparedMetadata, saveDir, FDRthreshold)
% Create violin plots of regression results
% INPUT:
% OUTPUT of:
% [pdRegressionResults,preparedInputTable,preparedMetadata] = performParkinsonAnalysis(dataPath, metadataPath, saveDir,'fluxRegressions.csv');

if nargin<5
    FDRthreshold = 0.11;
end

% Find regression results for metabolites with FDR < 0.1
if isstruct(regressionResults)
    fluxRes = regressionResults.Flux;
else
    fluxRes = regressionResults;
end
fluxRes(fluxRes.FDR>FDRthreshold,:)=[];

% Filter the preparedInputTable on metabolites with FDR<0.1.
% testedMets = preparedInputTable.Properties.VariableNames';
% idxToVis = matches(testedMets,fluxRes.Reaction);
% idxToVis(1) = true;
% fluxesToVis = preparedInputTable(:,idxToVis');
fluxesToVis = preparedInputTable(:,[{'ID'}, cellstr(fluxRes.Reaction)'] );


% Add disease status for visualising the differences between PD patients
% and controls.
metadataToVis = preparedMetadata(:,{'ID','Case_status'});
fluxesToVis = outerjoin(fluxesToVis,metadataToVis,'Type','left','MergeKeys',true);
fluxesToVis = movevars(fluxesToVis,'Case_status','After','ID');

% Process disease status groups
fluxesToVis.Case_status(matches(fluxesToVis.Case_status,'Control')) = {'NHC'};

% Translate VMH IDs to common metabolite names
metsToTest = string(fluxesToVis.Properties.VariableNames(3:end))';
metNames = string(renameVmhToMetName(metsToTest)); % Function that is also later used
metsToTest = [metsToTest metNames];
% metsToTest(matches(metsToTest,"DM_leu_L[bc]"),2) = "L-leucine";
% metsToTest(matches(metsToTest,"DM_leuleu[bc]"),2) = "Leucyleucine";
% metsToTest(matches(metsToTest,"DM_but[bc]"),2) = "Butyrate";
% metsToTest(matches(metsToTest,"DM_pnto_R[bc]"),2) = "Pantothenate";
% metsToTest(matches(metsToTest,"DM_nac[bc]"),2) = "Nicotinic acid";
% metsToTest(matches(metsToTest,"DM_ttdca[bc]"),2) = "Myristic acid";
size(metsToTest,1)

% Create figure
if size(metsToTest,1)==6
    
    fig = figure('Position',[-1663,144,1512,828]);
    t = tiledlayout(2,3,'TileSpacing','tight','Padding','loose');
    type = " fluxes in blood";
elseif size(metsToTest,1)==2

    fig = figure('Position',[258,293,1090,457]);
    t = tiledlayout(1,2,'TileSpacing','tight','Padding','loose');
    type = " associated microbial abundances";
end

for i=1:size(metsToTest,1)
    nexttile;
    plotViolin(fluxesToVis,metsToTest, fluxRes, i, type);
end


if ~matches(type," associated microbial abundances")
    ylabel(t,"Normalised log2 flux in mmol/day/person",'FontSize',18)
    xlabel(t,"Neurologically healthy controls (NHC) vs Parkinson's disease (PD) patients",'FontSize',18)

    % Save figure
    exportgraphics(fig,fullfile(saveDir,'fluxResultsPD.png'),'Resolution',300)
else
    ylabel(t,"Normalised log2 summed relative abundances",'FontSize',14)
    xlabel(t,"Neurologically healthy controls (NHC) vs Parkinson's disease (PD) patients",'FontSize',14)

    % Save figure
    exportgraphics(fig,fullfile(saveDir,'microbiomeSubsetsPD.png'),'Resolution',300)
end


end

function plt = plotViolin(fluxesToVis,metsToTest, fluxRes, i, type)

% Create grouped violin for metabolite i
plt = violinplot(fluxesToVis.(metsToTest(i,1)),fluxesToVis.Case_status);

% TODO: Update colours of violin plots

% Add axis labels and title
title(strcat(metsToTest(i,2), type),'FontWeight','normal')
%ylab = 'Normalised log2 flux in mmol/day/person';
ylabel('')
xlabel('')

% Add N samples to plot legend
ax=gca;
pdIdx = matches(ax.XTickLabel,'PD');
ax.XTickLabel(pdIdx==0) = append(ax.XTickLabel(pdIdx==0), ' (N=',cellstr(string( fluxRes.Control_N(i) )), ')'); % Add control sample size
ax.XTickLabel(pdIdx==1) = append(ax.XTickLabel(pdIdx==1), ' (N=',cellstr(string( fluxRes.PD_N(i) )), ')'); % Add PD sample size

% Format axis labels
set(gca,'FontSize',12)
set(gca,'TitleFontSizeMultiplier',1.5)
if matches(type," associated microbial abundances")
    set(gca,'TitleFontSizeMultiplier',1.1)
end

% Add regression p-values to plot
pVal = fluxRes.pValue(matches(fluxRes.Reaction,metsToTest(i,1)));

% Setup plot annotation with p-value
textInput = ['p = ',sprintf('%0.2E',pVal)];
xcoord = 1.5;
%ycoord = max(fluxesToVis.(metsToTest(i,1))) * 0.9;
yRange = max(ylim) - min(ylim);
ycoord = max(ylim) - yRange/25;
fSize = 14;

% Annotate plot with p-value
text(xcoord,ycoord,textInput,'HorizontalAlignment','center','FontName','Arial','FontSize',fSize,'Interpreter','none');

% Annotate plot with plot names
subPlotNames = 'abcdefgh';
titleProperties = get(gca,'Title');
text(0.5, titleProperties.Position(2),subPlotNames(i),'FontSize',12*1.5,'VerticalAlignment','bottom')
end
