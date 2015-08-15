% Compute from two matrices u, v such that X=U*V^T
% to matrices representing them when multiplied by each other 
% to get X in a vectorized form 
%
% Input:
% U - matrix size n1Xr 
% V - matrix size n2Xr
%
% Output:
% U_mat - matrix of size (n1*n2)X(r*n2) - such that U_mat*vec(V) = vec(X)
% V_mat - matrix of size (n1*n2)X(r*n2) - such that V_mat*vec(U) = vec(X)
% x_vec - vector obtained by stacking the coluns of X=U*V^T
%
function [u_mat, v_mat, x_vec] = bilinear_to_general_affine(U, V)
[n1, r] = size(U);
[n2, r] = size(V);
u_mat = kron(eye(n2), U);
v_mat = kron(V, eye(n1));
X = U*V';
x_vec = mat2vec(X);

