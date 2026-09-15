function A =construct_linearApproximation(m,n,ii,N,X,Y,x_L,x_U,gamma)
    D = zeros(N,N);
    Delta_X = zeros(m,N);
    Delta_Y = zeros(n,N);
    for jj =1:N
        if jj~=ii
            value = 0;
            for kk =1:m
                value = value + ((X(kk,jj)-X(kk,ii))/(x_U(kk)-x_L(kk)))^2;
            end
            D(jj,jj) = (1/value)^gamma;
        else
            D(jj,jj) = 0;
        end
        Delta_X(:,jj) = X(:,jj)-X(:,ii);
        Delta_Y(:,jj) = Y(:,jj)-Y(:,ii);
    end
    A = (Delta_Y*D)*pinv(Delta_X*D);
end