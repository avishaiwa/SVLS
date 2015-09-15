% Solve the following least-squares problem: 
% x_hat = argmin_x || [A (*) B] x - b ||_2^2 
% where (*) denotes the matrix Kronecker product 
% 
% Input: 
% A1 - first matrix in Kronecker product for design matrix
% A2 - second matrix in Kronecker product for design matrix
% b - vector of measurements
% alg_str - algorithm to use 
% 
% Ouptut: 
% x - least-squares solution 
% 
function x = Kronecker_LS(A1, A2, b, alg_str)

if(~exist('alg_str', 'var') || isempty(alg_str))
    alg_str  = 'kronecker_psuedo_inv';
end


switch lower(alg_str)
    case {'brute_force', 'ls', 'least_squares'} % don't use kronecker product. Solve standard least-squares problem. 
        x = kron(A1, A2) \ b; 
    case 'kronecker_psuedo_inv' % Compute psuedo-inverses for A1, A2
        x = mat2vec( pinv(A2) * vec2mat(b', size(A2,1))' * (pinv(A1))' ); 

    case 'kronecjer_qr_decomposition' % Use full method based on QR-decomposition in paper by [Fausett & Fulton 1994]     
    % T.B.D. 
        
end
