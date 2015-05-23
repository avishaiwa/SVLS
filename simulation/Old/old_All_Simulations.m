% Run simulations for low-rank matrix recovery
%
% Input:
% n - size of matrix X to be sampled (assume square: nXn)
% r - rank of matrix X.
% k_iter - number of sampled matrix for each k number of measurements
% alg_str - vector of recovery algorithms
% maximum_iter - maximum iterations in gradient descent
% tol - stopping rule for the gradient descent
% index -  number of columns and rows (or measurements) that are sampled
% X_type - Type of unknown matrix X: low_rank or symmetric_low_rank
% measurement_type - uniform, columns_and_rows or normal
%
% Output:
% losses, times and iterations are 3-dimentional arrays that store the RMSE,
% the time and number of iteration for each run, respectivelly :
% three - index of the algorithm in alg_str
% two - number of measurement (p or k)
% one - number of iterations (i)
% For example: losses(i,p,a) is the rmse obtained by running algorithm a
% with p measurements on the i-th run
%
function [ losses, times, iterations ] = All_Simulations( ...
    n,r,noise,alg_str,k_iter,maximum_iter,tol,index ,X_type,measurement_type)


if(~exist('tol', 'var') || isempty(tol))
    tol = 0.0001; % set default tolerance
end
if(~exist('maximum_iter', 'var') || isempty(maximum_iter))
    maximum_iter = 100; % set default # of iterations
end
if(~exist('k_iter', 'var') || isempty(maximum_iter))
    k_iter= 1; % set default # of iterations
end
if(~exist('X_type', 'var') || isempty(maximum_iter))
    X_type = 'low_rank';% set default matrix type
end
if(~exist('measurement_type', 'var') || isempty(maximum_iter))
    measurement_type = 'columns_and_rows'; % set default measurements type
end


num_algs = length(alg_str);
losses = zeros(k_iter, length(index), num_algs); 
iterations = losses; times = losses; 
ind=0;

for p=index% loop on # of measurements
    ind = ind+1    
    for i=1:k_iter % loop on iterations per measurement
        [m, am, Br, Bc, Ac, Ar, map, tmap] = ...
            samp_matrix(n, r, p, noise, X_type, measurement_type);
        
        for t = 1:num_algs
            [~, time, loss, iteration] = AffineMatrixRecovery( ...
                m, am, alg_str{t}, p,noise, maximum_iter, tol, map, tmap, Br, Bc, Ac, Ar);
            
            losses(i,ind,t)=loss;
            times(i,ind,t)=time;
            iterations(i,ind,t)=iteration;
                        
        end
    end
end


