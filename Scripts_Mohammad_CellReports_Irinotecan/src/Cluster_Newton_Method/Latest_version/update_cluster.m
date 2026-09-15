function [X_k1,Y_k1,r_k1,lamda_k1] = update_cluster(m,n,N,f,X_k,X_k1,Y_k,Y_star,r_k,lamda_k,lamda_max)
    Y_k1 = zeros(n,N);
    r_k1 = zeros(1,N);
    lamda_k1 = zeros(1,N);
    parfor ll =1:N
        dummy = f(X_k1(:,ll));
        if ((lamda_k(ll) <=  lamda_max) & all(~isnan(dummy)) & all((dummy)>0) & all((X_k1(:,ll))>0))
            Y_k1(:,ll) = dummy;
            %r_k1(ll) = norm(Y_k1(:,ll)-Y_star);
            r_k1(ll)= norm(diag(Y_k1(:,ll)-Y_star)*[1./Y_star]);
            if r_k1(ll) > r_k(ll)
                X_k1(:,ll) = X_k(:,ll);
                Y_k1(:,ll) = Y_k(:,ll);
                r_k1(ll) = r_k(ll);
                lamda_k1(ll) = 10* lamda_k(ll);
            else
                lamda_k1(ll) = lamda_k(ll)/10;
            end
        else
            X_k1(:,ll) = X_k(:,ll);
            Y_k1(:,ll) = Y_k(:,ll);
            r_k1(ll) = r_k(ll);
            lamda_k1(ll) = 10* lamda_k(ll);
        end
    end
end 