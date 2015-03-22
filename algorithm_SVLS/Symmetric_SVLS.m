%try to reconstract symmetric matrix m from rank r and size nXn with SVLS algorithm.
% Input: 
% B - matrix of row measuremts of size k * n
% A - matrix of row measuremts coefficients of size k * n
% r - rank of unknown matrix 
% 
% Output: 
% X - the estimate of the unknown symmetric matrix M where we minimize  ||M*ac-Bc||^2 +  ||ar*M-Br||^2 
%
function [ X, est_r] = Symmetric_SVLS( B, A, r )

if (~exist('r', 'var') || isempty(r))
    fprintf(1,'Rank not specified. Trying to guess ...\n');
    s = svd(B);
    [~, r] = max(s(1:end-1)./s(2:end)); 
end
est_r = r

addpath('PROPACKmod'); % for top-svd algorithm 
[u, ~, ~] = lansvd(B', r,'L','OPTIONS'); %find the r largest vectors and singular values.
X = (u*(A*u\B))';
