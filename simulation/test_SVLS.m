% Sample matrix and run SVLS algorithm 
n=100; r=6; noise=0.95; X_type='low_rank'; measurement_type = 'gaussian_columns_and_rows';
alg_str = {'svls_alternating_ls'}; p_iter=5; max_iter=100; tol=10^-5; k_vec=24;

[ losses, times, iterations] = All_Simulations( ...
    n,r,noise,alg_str,p_iter,max_iter,tol, k_vec ,X_type, measurement_type);

fprintf('Sampled matrix: n=%d, rank=%d, noise: sigma=%f\n', n, r, noise)
fprintf('Runnning SVLS with k=%d row-column measurements (%d iterations)\n', k_vec, p_iter)
fprintf('Average RMSE Error: %f\n', mean(losses))
