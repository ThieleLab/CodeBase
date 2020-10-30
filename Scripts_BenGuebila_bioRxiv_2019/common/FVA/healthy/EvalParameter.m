function value = EvalParameter(P, t, y)

    if isnumeric(P)
        value = P;
    else
        value = P(t,y);
    end
