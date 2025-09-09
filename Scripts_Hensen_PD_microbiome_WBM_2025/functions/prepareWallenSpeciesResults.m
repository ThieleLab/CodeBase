function speciesPDdata = prepareWallenSpeciesResults(parkinsonMicrobes, savePath)
% Process gut microbial species shifts in PD patients associated in Wallen et al. (2022)

% INPUT:
% path to annotated species data

% OUTPUT
% path to saved processed species data with matching species names 

% Load dataset
speciesPDdata = readtable(parkinsonMicrobes,'VariableNamingRule','preserve');

% Remove p-value data
speciesPDdata(:,'MaAsLin2 P')=[];
speciesPDdata(:,'ANCOM-BC P')=[];

% Rename species columns and move
speciesPDdata = renamevars(speciesPDdata,'Species','Species_Wallen');
speciesPDdata = renamevars(speciesPDdata,'Species_AGORA2APOLLOtaxonomy','Species');
speciesPDdata = renamevars(speciesPDdata,'Famiy','Family');

speciesPDdata = movevars(speciesPDdata,'Species','After','Species_Wallen');

% Get significantly changed PD species
speciesPDdata.Significant = speciesPDdata.("MaAsLin2 FDR")<0.1 & speciesPDdata.("ANCOM-BC FDR")<0.1;

% Get direction of change
speciesPDdata.Direction = repmat("", height(speciesPDdata), 1);
speciesPDdata.Direction(speciesPDdata.Significant == 1 & speciesPDdata.("MaAsLin2 Beta")>0 & speciesPDdata.("ANCOM-BC Beta")>0) = 'Higher in PD';
speciesPDdata.Direction(speciesPDdata.Significant == 1 & speciesPDdata.("MaAsLin2 Beta")<0 & speciesPDdata.("ANCOM-BC Beta")<0) = 'Lower in PD';
speciesPDdata.Direction(speciesPDdata.Significant == 0) = 'Unchanged';

% Save processed microbial species results
writetable(speciesPDdata,savePath)
end