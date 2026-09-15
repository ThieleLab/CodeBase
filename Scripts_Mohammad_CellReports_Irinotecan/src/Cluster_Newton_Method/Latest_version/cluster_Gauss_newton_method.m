function [X,Data] = cluster_Gauss_newton_method(m,n,x_L,x_U,N,f,Y_star,lamda_init,lamda_max,gamma,k_max)
    r = [ones(1,N)];
    for ii=1:N
        lamda(ii) = lamda_init;
    end 
    [X,Y,r] = create_initial_cluster(m,n,x_L,x_U,N,f,Y_star);    
    for k= 0: k_max-1
        disp(['Iteration Number: ',num2str(k)])
        Data(k+1,:)={k,X,Y,r,lamda,[]};
        for ii =1:N
            if lamda(ii)<=lamda_max
                A(:,:,ii) = construct_linearApproximation(m,n,ii,N,X,Y,x_L,x_U,gamma);
                X_k1(:,ii) = X(:,ii)+(transpose(A(:,:,ii))*A(:,:,ii)+lamda(ii)*eye(m))^-1*transpose(A(:,:,ii))*(Y_star-Y(:,ii));
            else 
            end
        end
        Data(k+1,end)={A};
        [X,Y,r,lamda] = update_cluster(m,n,N,f,X,X_k1,Y,Y_star,r,lamda,lamda_max);
    end
end 