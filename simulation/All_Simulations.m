% Run simulations for low-rank matrix recovery
% Simulate matrices and measurements, recover them from measurements using
% an algorithm of choice, and record algorithm's performance
%
% Input:
% n - size of matrix X to be sampled (assume square: nXn)
% r - rank of matrix X.
% num_sampled_matrices - number of sampled matrices for each k number of measurements
% alg_str - vector of recovery algorithms
% maximum_iter - maximum iterations in gradient descent
% tol - stopping rule for the gradient descent
% k_vec -  number of columns and rows (or measurements) that are sampled
% X_type - Type of unknown matrix X: low_rank or symmetric_low_rank
% measurement_type - uniform, columns_and_rows or normal
% known_rank - 0: assume rank is unknown and estimate it (default). 1: assume rank is known
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
    n, r, noise, alg_str, num_sampled_matrices, maximum_iter, tol, k_vec , ...
    X_type, measurement_type, known_rank)


if(~exist('tol', 'var') || isempty(tol))
    tol = 0.0001; % set default tolerance
end
if(~exist('maximum_iter', 'var') || isempty(maximum_iter))
    maximum_iter = 100; % set default # of iterations
end
if(~exist('num_sampled_matrices', 'var') || isempty(num_sampled_matrices))
    num_sampled_matrices= 1; % set default # of iterations
end
if(~exist('X_type', 'var') || isempty(X_type))
    X_type = 'low_rank';% set default matrix type
end
if(~exist('measurement_type', 'var') || isempty(measurement_type))
    measurement_type = 'columns_and_rows'; % set default measurements type
end
if(~exist('known_rank', 'var') || isempty(known_rank))
    known_rank = 0; % default is estimating rank from data 
end
if(~known_rank) % set unknown rank 
    r = []; 
end


num_algs = length(alg_str);
losses = zeros(num_sampled_matrices, length(k_vec), num_algs); 
iterations = losses; times = losses; 

for j=1:length(k_vec) % loop on # of measurements
    simulate_k = k_vec(j) % ind = ind+1    
    for i=1:num_sampled_matrices % loop on iterations per measurement
        sprintf('Run k=%ld, iter=%ld out of %ld\n', k_vec(j), i, num_sampled_matrices) 
        [X, measure] = ...
            samp_matrix(n, r, k_vec(j), noise, X_type, measurement_type);
        
        for t = 1:num_algs
            [~, times(i,j,t), losses(i,j,t), iterations(i,j,t)] = AffineMatrixRecovery( ...
                X,  alg_str{t}, k_vec(j), noise, maximum_iter, tol, measure, r);                        
        end
    end
end


