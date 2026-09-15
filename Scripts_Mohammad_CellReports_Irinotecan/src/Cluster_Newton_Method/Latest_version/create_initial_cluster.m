function [X,Y,r] = create_initial_cluster(m,n,x_L,x_U,N,f,Ystar)

X = [NaN*ones(m,N)]; % allocate memory for mXN matrix
Y = [NaN*ones(n,N)]; % allocate memory for nXN matrix
r = [ones(1,N)]; % 1XN vector

parfor ll =1:N
    while all(isnan(Y(:,ll)))==1
        X(:,ll) = unifrnd(x_L,x_U,m,1); % Uniformlt generate cluster 
        dummy = f(X(:,ll));
        if (all(~isnan(dummy))==1 & all((X(:,ll))>0))
            Y(:,ll) = dummy;
        else
            Y(:,ll) = [NaN*ones(n,1)];
        end
    end
%   Sum of sqaured residual 
    %r(:,ll)=norm(Y(:,ll)-Ystar);
    r(:,ll)=norm(diag(Y(:,ll)-Ystar)*[1./Ystar]);
end

end