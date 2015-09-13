% Sample matrix and run SVLS algorithm 
tol = 10^(-15); 
n=30; r=6; noise=0.0; X_type='low_rank'; measurement_type= 'gaussian_columns_and_rows'; % 'columns_and_rows';
alg_str = {'svls'}; p_iter=5; max_iter=100; k=5; % tol=10^-5; k=5;

[ X, measure ] = samp_matrix( ...
    n,r,k,noise, X_type, measurement_type);

[AA, bb] = row_column_to_general_affine(measure.Ar, measure.Ac, measure.Br, measure.Bc);

figure; imagesc(AA); colorbar;

sum(AA)

% Test conversion
find(measure.Ar*X - measure.Br)
find(X*measure.Ac - measure.Bc)
vec_b = AA*mat2vec(X);
x_vec = mat2vec(X); 
bad_inds = find(abs(AA*(mat2vec(X)) - bb')>tol) % should be empty !!! 

% Create 'artificial 
[U, Lambda, V] = svd(X); U = U*Lambda; U=U(:,1:r); V=V(:,1:r); 
max(max(abs(U*V' - X)))
[U_mat, V_mat, x_vec2] = bilinear_to_general_affine(U, V);
find( abs(U_mat * mat2vec(V') - x_vec2) > 0.00000001 )
find( abs(V_mat * mat2vec(U) - x_vec2) > 0.00000001 )

