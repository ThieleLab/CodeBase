
% Very time-consuming. You may need to restart multiple times as MATLAB
% sometimes runs out of memory during the creation of SBML files.
dInfo = dir(refinedFolder);
modelList={dInfo.name};
modelList=modelList';
modelList(~contains(modelList(:,1),'.mat'),:)=[];

mkdir([pwd filesep 'SBML_Files'])
% remove already created
dInfo = dir([pwd filesep 'SBML_Files']);
sbmlList={dInfo.name};
sbmlList=sbmlList';
sbmlList(~contains(sbmlList(:,1),'.xml'),:)=[];
sbmlList = strrep(sbmlList,'.xml','.mat');
modelList = setdiff(modelList,sbmlList); 

for i=1:length(modelList)
    model=readCbModel([refinedFolder filesep modelList{i}]);
    writeCbModel(model,'format','sbml','fileName',[pwd filesep 'SBML_files' filesep strrep(modelList{i},'.mat','')])
end