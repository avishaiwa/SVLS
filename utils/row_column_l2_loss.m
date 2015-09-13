% L2 loss function for decomposition
%
% Input:
% W - First matrix in decomposition of X_hat (X_hat = W*S)
% S - second matrix in decomposition of X_hat
% Br = Ar*X + Normal(0,noise)
% Bc = X*Ac + Normal(0,noise)
% Ar = matrix size pXn - row measurements
% Ac = matrix size nXp - column measurements
%
% Output:
% L2_loss - Loss F(W,S) = ||Ar*W*S'-Br||_F^2 + ||W*S'*Ac-Bc||_F^2|
%
function L2_loss = row_column_l2_loss(W, S, Br, Bc, Ar, Ac)

L2_loss =0.5*norm(Ar*W*S-Br,'fro')^2+0.5*norm(W*S*Ac-Bc,'fro')^2;
