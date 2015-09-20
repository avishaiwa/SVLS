% Compute the product: y = [A (*) B] x 
% where (*) denotes the matrix Kronecker product 
% 
% Input: 
% A - first matrix in Kronecker product for multiplier matrix
% A - second matrix in Kronecker product for multiplier matrix
% x - vector to multiply by 
% alg_str - algorithm to use 
% 
% Ouptut: 
% y - product: y = [A (*) B] x 
% 
function y = Kronecker_matrix_vector_product(A, B, x, alg_str)

if(~exist('alg_str', 'var') || isempty(alg_str))
    alg_str  = 'kronecker_prod';
end

switch lower(alg_str)
    case {'brute_force'} % don't use kronecker product. Solve standard least-squares problem. 
        y = kron(A, B) * x; 
    case 'kronecker_prod' % Compute psuedo-inverses for A1, A2
        y = mat2vec( B * vec2mat(x', size(B,2))' * A' ); 
end


