function model=augmentHarvey(model,nCores);
% Adds metasusbystems 
%      organ vector 
%      gene rxn matrix
%      metasubsystem vector
%      recon reaction vector
%      recon metabolite vector
% Marouen BEN GUEBILA 01/2016
    n = length(model.rxns);
    %%
    %fix subsystem length
    if length(model.subSystems) < length(model.rxns)
        model.subSystems(end-(length(model.rxns)-length(model.subSystems))+1:end)='';
    end
    %%
    %add gene rxn matrix to harvey
    part = floor(n/nCores);
    matSave = cell(nCores,1);
    parfor j=1:nCores
        harvey=model;
        for i=((j-1)*part)+1:j*part
            harvey = changeGeneAssociation(harvey,harvey.rxns{i},harvey.grRules{i});
        end
        matSave{j} = harvey.rxnGeneMat;
    end
    catA=cat(nCores,matSave{:});
    model.rxnGeneMat = sum(catA,nCores);
    if mod(n,nCores)
        for i = nCores*part:n
            model = changeGeneAssociation(model,model.rxns{i},model.grRules{i});
        end
    end
    %%
    %add organ vector to harvey
    organName = cellfun(@(x) x(1:strfind(x,'_')-1),model.rxns,...
        'UniformOutput',false);
    %Manually delete compartments
    comps= {'BileDuct' 'Diet' 'EX' 'Excretion' 'Whole' 'sink' 'BBB' 'GI' 'LI' 'SI'};
    organName = unique(organName(cell2mat(cellfun(@(x) ~ismember(x,comps),organName,...
        'UniformOutput',false))));
    model.organs = cell(n,1);

    for organ = organName'
        l = length(find(cellfun(@(x) ~isempty(strfind(x,[organ{1} '_'])),model.rxns)));
        model.organs(find(cellfun(@(x) ~isempty(strfind(x,[organ{1} '_'])),model.rxns))) = repmat(organ,[1 l]); 
    end
    
    lNa=length(find(cellfun(@(x) isempty(x),model.organs)));
    model.organs(find(cellfun(@(x) isempty(x),model.organs))) = repmat({'NA'},[1 lNa]);
    %%
    %add metasubsystem vector to harvey
    model.metaSubSystems = cell(n-1,1);
    metaSub = readtable('metaSubsystems.csv','ReadVariableNames',false);

    for i = 1:size(metaSub(2:end,1),1)
       subSystemStr = table2cell(metaSub(i,1));
       subSystemStr = subSystemStr{:};
       MetaSubStr = table2cell(metaSub(i,2));
       l=length(find(cellfun(@(x) ~isempty(strfind(x,subSystemStr(2:end-1))),model.subSystems)));
       model.metaSubSystems(find(cellfun(@(x) ~isempty(strfind(x,subSystemStr(2:end-1))),model.subSystems))) = repmat(MetaSubStr,[1 l]);
    end
    lNa=length(find(cellfun(@(x) isempty(x),model.metaSubSystems)));
    model.metaSubSystems(find(cellfun(@(x) isempty(x),model.metaSubSystems))) = repmat({'NA'},[1 lNa]);
    model.metaSubSystems(end+1)={'NA'};
    %%
    %recon reaction vector
    model.uniqueRxn = cellfun(@(x) char(regexp(x, '(?<=_)\w+((?=[))?','match')),model.rxns,'UniformOutput',false);
    %%
    %compaitibility with functions using A or S
    model.Smat = model.S;
    model.S = model.A;
    %%
    %recon metabolite vector
    a=cellfun(@(x) char(regexp(x, '\w+((?=[))','match')),model.mets,'UniformOutput',false);
    a=a(1:size(model.Smat,1));
    for i=unique(model.organs)'
        if isequal(i{1},'NA')
            continue
        end
        i
        b=find(cellfun(@(x) ~isempty(strfind(x,[i{1} '_'])),a));
        c=cellfun(@(x) char(regexp(x, '(?<=_)\w*','match')),a(b),'UniformOutput',false);
        a(b)=c;
    end
    model.uniqueMet=a;
end