function newVec = vecBuilder(oldNumVec,harvey,stochvec);
    [oldNumVec,indexSort] = sort(oldNumVec);
    %sort stochiometry vector by order of oldNumVec
    stochvec = stochvec(indexSort);
    %initilaize newvector
    newVec = [zeros(1,oldNumVec(1)-1) stochvec(1)];
    i = 1;
    for i = 2:length(oldNumVec)
        vec = [zeros(1,(oldNumVec(i)-1)-(oldNumVec(i-1)+1)+1) stochvec(i)];
        newVec = [newVec vec];
    end
    if isequal(i,[])
        i = 1;
    end
    if i == length(oldNumVec)
        vec = zeros(1,length(harvey.rxns)-(oldNumVec(i)+1)+1);
        newVec = [newVec vec];
    end
end