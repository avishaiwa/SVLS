% Relative Root-Mean-Squared-Error for two matrices
% Input:
% m-the matrix that sampeled
%
% Output:
% out - RR - || X1-X2||_F / ||X2||_F
function RR = RRMSE(X1,X2)
RR = norm(X1-X2,'fro')/norm(X2,'fro');

