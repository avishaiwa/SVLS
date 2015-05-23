% Test reconstruction performance of basis algorithm using simulations
AssignGeneralConstants;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters: 
n = 50; % set matrix dimension
r = 10; % set rank 
sigma_noise =0.2; % set noise level 
alg_vec = {'basis_gradient1', 'OptSpaceNoTrim'};
p_iter = 10; % number of simulations to perform at each #of measurements
gradient_stepsize = 5*10^-4;
max_iter = 100; % maximum iterations in gradient descent
tol=10^(-4); % tolerane for gradient descent algorithm
min_k = 15; % set minimal # of row&column measurements. 
max_k = 20; % set number of row and column measurements. 
            % The number of scalar measurements is k*2*n
X_type =  'low_rank';
measurement_type = 'columns_and_rows';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulations to get errors and times 
[ losses, times, iterations ] = All_Simulations( ...
    n,r,sigma_noise,alg_vec,p_iter,max_iter,tol,min_k,max_k ,X_type,measurement_type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results and save to files: 
% Plot reconstuction error as function of # samples (need to loop on algorithms)
figure; boxplot(losses(:,min_k:max_k,2), 'labels', min_k:max_k); 
xlabel('Num. row-column measurements'); ylabel('RMSE error'); 


% Replace boxplot with error-bar
figure; hold on; epsilon=0.1; 
for a=1:length(alg_vec)
    errorbar((min_k:max_k) + epsilon*(a-1), ...
        mean(losses(:,min_k:max_k,a)), std(losses(:,min_k:max_k,a)), color_vec(a)); 
end
legend(alg_vec); legend('boxoff'); 
xlabel('Num. row-column measurements'); ylabel('RMSE error'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% New test: keep r and k fixed, and increase n. Does RMSE go to zero? (it shouldn't!)

n_vec = ceil(logspace(1.5,4,30));
p_iter = 5; % number of iterations per n 
r = 5; k=30; % set 
sigma_noise = 1; % set noise level 
RMSE_mat = zeros(p_iter, length(n_vec)); 
RMSE_vec = zeros(1,length(n_vec)); 

for i=1:length(n_vec)
    run_i = i 
    run_n = n_vec(i) 
    [ losses, times, iterations ] = All_Simulations( ...
        n_vec(i),r,sigma_noise, {'basis_gradient1'}, p_iter,0,tol, ...
        k,k ,X_type,'columns_and_rows_normal');
    RMSE_vec(i) = mean(losses(:,end)); 
    RMSE_mat(:,i) = losses(:,end); 
end    

figure; loglog(n_vec, RMSE_vec, '*'); 
AB = polyfit(log(n_vec), log(RMSE_vec), 1); 
hold on; plot(n_vec, exp(AB(1)*log(n_vec) + AB(2)), 'r'); 
figure; errorbar(mean(RMSE_mat), std(RMSE_mat)); 
figure; boxplot(RMSE_mat, 'labels', n_vec);
% Fit on loglog scale 

% TODO: 
% 1. Make plots in a function 
% 2. Make all input/output in a few structs 




