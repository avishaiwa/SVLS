% Reconstract matrix m of rank r and size nXn with SVLS algorithm.
% Input:
% Br - matrix of row measurements (Br = Ar*m + noise*N(0,1))
% Bc - matrix of column measurements (Bc = m*Ac + noise*N(0,1))
% Ar - row design matrix
% Ac - column design matrix
% r - rank of unknown matrix (optional) 
%
% Output:
% X - estimate of M where we minimize  ||m*ac-Bc||^2 +  ||ar*M-Br||^2
% est_r - estimated rank (when rank r is unknown)
%
function [X, est_r] = SVLS( Br,Bc,Ar,Ac,r)
if (~exist('r', 'var') || isempty(r))
    %  fprintf(1,'Rank not specified. Trying to guess ...\n');
    s = svd(Bc);
    [~, r1] = max(s(1:end-1)./s(2:end)); 
     s = svd(Br);
    [~, r2] = max(s(1:end-1)./s(2:end)); 
    r = max(r1,r2);
end
est_r = r;

[uC, ~, ~] = lansvd(Bc,r,'L','OPTIONS'); %find the r largest vectors and singulr values of Bc.
[uR, ~, ~] = lansvd(Br',r,'L','OPTIONS'); %find the r largest vectors and singulr values of Br.
Xr=(uR*(Ac'*uR\Bc'))';% find Xr
Xc=uC*(Ar*uC\Br); % find Xc

X = (Xr+Xc)/2;
%[u ,d, v]  = lansvd(X,r,'L','OPTIONS');
%X = u(:,1:r)*d(1:r,1:r)*v(:,1:r)';
end

