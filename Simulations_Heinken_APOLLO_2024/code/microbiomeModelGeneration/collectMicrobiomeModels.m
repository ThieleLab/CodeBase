% collect all created models in one folder

clear all
rootDir = pwd;

modPath = [rootDir filesep 'data' filesep 'GutMicrobiomes' filesep 'MicrobiomeModels'];
mkdir(modPath)
for i=1:14
    resPath=[rootDir filesep 'data' filesep 'GutMicrobiomes' filesep 'MicrobiomeModels_' num2str(i)];
    dInfo = dir(resPath);
    modelList={dInfo.name};
    modelList=modelList';
    modelList(~contains(modelList(:,1),'microbiota_model'),:)=[];
    for j=1:length(modelList)
        copyfile([resPath filesep modelList{j}],[modPath filesep modelList{j}])
    end
    % delete the previous folder
    rmdir(resPath, 's')
end
