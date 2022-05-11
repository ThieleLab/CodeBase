
% get 10 random AGORA2 reconstructions to test in MEMOTE.

infoFile = readInputTableForPipeline('AGORA2_infoFile.xlsx');
infoFile(1,:) = [];

agoraPath = [rootDir filesep 'Current_Version_AGORA2' filesep 'Output_Models'];

models = infoFile(randperm(length(infoFile),10));

mkdir('ReconsForMemote')
cd('ReconsForMemote')

for i=1:length(models)
 model = readCbModel([agoraPath filesep models{i} '.mat']);
 writeCbModel(model,'format','sbml','fileName',models{i})
end