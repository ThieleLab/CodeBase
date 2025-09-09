function dietTable = collectDietInfo(mWbms)
% input 
% mWbms: Path to folder with prepare WBMs.

% output
% dietTable: Table WBM diet setup

mdls = what(mWbms).mat; % Get model names
mdlRxnPar = load(fullfile(mWbms,mdls{1}),'rxns','lb','ub'); % Read host-microbiome WBM reactions and the set reaction bounds.
mdlRxnPar = struct2table(mdlRxnPar);
% Get the dietary reactions with negative lower bounds (metabolites that
% can be consumed.
dietContent = mdlRxnPar(contains(mdlRxnPar.rxns,'Diet_') & mdlRxnPar.lb<0 , :);

%
mets = erase(dietContent.rxns,{'Diet_EX_','[d]'}); % Extract VMH metabolites
%mets(matches(mets,'adpcbl')) = {'adocbl'}; % Fix wrong metabolite name in VMH database
DB = loadVMHDatabase().metabolites; % Load VMH database

[~,~,ib]=intersect(mets,DB(:,1),'stable'); % Find metabolite names
metMapping = array2table(DB(ib,[1 2]),'VariableNames',{'rxns','Metabolite'});
% Unfix adocbl for consistency 
%metMapping.rxns(matches(metMapping.rxns,'adocbl')) = {'adpcbl'};
metMapping.rxns = append('Diet_EX_',metMapping.rxns,'[d]');

% Add metabolite names to dietContent information
dietTable = outerjoin(dietContent,metMapping,'Keys','rxns','MergeKeys',true);

% Manually apply fix 
dietTable.Metabolite(matches(dietTable.rxns,'Diet_EX_adpcbl[d]')) = {'Adenosylcobalamin'};

% Prepare table for supplementary materials
dietTable = movevars(dietTable,'Metabolite','After','rxns');
dietTable.Properties.VariableNames = {'Diet reaction',...
    'Associated metabolite',...
    'Lower flux bound in mmol/day/person',...
    'Upper flux bound in mmol/day/person'};
end