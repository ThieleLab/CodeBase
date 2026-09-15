tic

m = 7; % Number of parameter
n = 3; % Number of objective function 
x_L = [0.006;0.006;0.006;0.006;0.006;0.006;0.006]; % Initial estimate lower bound
x_U =[0.6;0.6;0.6;6;0.6;6;0.6]; % Initial estimate Upper bound
N = 250; % Cluster size
f = @(X) (IRI_param_1(X)); % Irinotecan parameter defined for 1st patient
Y_star =  [23.2;45.3;1.53]; % Objective value for 1st patient
lamda_init = 0.01; % Minimum step size
lamda_max = 1e10; % Maximum step size
gamma = 1;
k_max = 250; % Number of iteration
[X,Data] = cluster_Gauss_newton_method(m,n,x_L,x_U,N,f,Y_star,lamda_init,lamda_max,gamma,k_max); % solve using cluster netwon method

toc