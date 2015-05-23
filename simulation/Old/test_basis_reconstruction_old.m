% Test reconstruction performance of basis algorithm 
n = 1000; % set matrix dimension
r = 10; % set rank 
k =20; % set number of row and column measurements 
sigma_noise =2; % set noise level 
gradient_stepsize = 5*10^-4;


[X, Br, Bc, Ar, Ac] = ColandRow_basis(n, r, k, sigma_noise,'symmetric_low_rank');
%    'symmetric_low_rank'); % simulate data for the symmetric case 

X_basis = basis(Br, Bc, Ar, Ac, r); % Solve using the basis algorithm (just svd part, no gradient) 
X_hat = gradientDescent1(X_basis, Br, Bc, Ar, Ac, r, [], gradient_stepsize); % Solve using the basis algorithm with gradient descent
X_basis_Symmetric = Symmetric_basis(Br, Ar, r);
B=0.5*(X_basis_Symmetric+X_basis_Symmetric'); % symmetrize matrix
X_basis_Symmetric2 = 0.5*(B+real(sqrtm(B'*B)));
X_hat_Symmetric = gradientDescent1(X_basis_Symmetric, Br, Bc, Ar, Ac, r, [], gradient_stepsize); % Solve using the basis algorithm with gradient descent
X_hat_Symmetric2 = gradientDescent1(X_basis_Symmetric2, Br, Bc, Ar, Ac, r, [], gradient_stepsize); % Solve using the basis algorithm with gradient descent

X_hat_Symmetric3 = Symmetric_gradientDescent( X_basis_Symmetric2, Br, Ar, r, [], gradient_stepsize); % run symmetric gradient descent

RMSE_basis = norm(X_basis-X,'fro')/norm(X,'fro') % error before gradient descent stage 
RMSE_basis_Symmetric = norm(X_basis_Symmetric-X,'fro')/norm(X,'fro')  % error after symmetric basis algorithm  
RMSE_basis_Symmetric2 = norm(X_basis_Symmetric2-X,'fro')/norm(X,'fro')   % error after symmetric basis algorithm  
RMSE_hat = norm(X_hat-X,'fro')/norm(X,'fro') % error after gradient descent stage 
RMSE_hat_Symmetric = norm(X_hat_Symmetric-X,'fro')/norm(X,'fro')   % error after gradient descent stage 
RMSE_hat_Symmetric2 =  norm(X_hat_Symmetric2-X,'fro')/norm(X,'fro')   % error after gradient descent stage 
RMSE_hat_Symmetric3 =norm(X_hat_Symmetric3-X,'fro')/norm(X,'fro') % error after symmetric gradient descent stage 


 [ am, X1] = opt_measurment(n, r, 2*n*k-k^2, sigma_noise);%, 'symmetric_low_rank');

[U_opt S_opt V_opt]  = optSpace(am,r,100,5*10^-2);
X1_opt = U_opt*S_opt*V_opt';
X1_opt_symetric = (X1_opt' +X1_opt )/2;
[U S V] = lansvd(X1_opt_symetric,r,'L','OPTIONS');
X1_opt_symetric2=U*S*V';

RMSE_opt = norm(X1_opt-X1,'fro')/norm(X1,'fro')   % error after gradient descent stage 
RMSE_opt_Symmetric2 =  norm(X1_opt_symetric-X1,'fro')/norm(X1,'fro')   % error after gradient descent stage 
RMSE_opt_Symmetric3 =norm(X1_opt_symetric2-X1,'fro')/norm(X1,'fro') % error after symmetric gradient descent stage 




