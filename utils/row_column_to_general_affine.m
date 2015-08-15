% Create a big matrix and measurement vector from measurement matrices in
% row-column model
%
% Input:
% Ar - matrix size kRXn1 - row measurements
% Ac - matrix size n1XkC - column measurements
% Br - matrix of size kRXn2 (equal to Ar*X + noise)
% Bc - matrix of size n2XkC (equal to X*Ac + noise)
%
% Output:
% X - estimation of m without any noise
% dis - matrix of distances (errors, divided into two terms) as function of iteration
%
function [AA, bb] = row_column_to_general_affine(Ar, Ac, Br, Bc)
[kR, n1] = size(Ar); % get first dimensions
[n2, kC] = size(Ac);

% Create big matrices with Kronecker products
AA = [kron(eye(n2), Ar)' kron(Ac', eye(n1))']'; % make big measurements matrix

%bb = [mat2vec(Br)' mat2vec(Bc)']; % vectorize measurements
bb = [Br(:)' Bc(:)']; 
end


