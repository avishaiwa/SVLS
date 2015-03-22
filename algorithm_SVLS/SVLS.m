% Reconstract matrix m of rank r and size nXn with SVLS algorithm.
% Input:
% Br - matrix of row measurements (Br = Ar*m + noise*N(0,1))
% Bc - matrix of column measurements (Bc = m*Ac + noise*N(0,1))
% Ar - row design matrix
% Ac - column design matrix
% r - rank of unknown matrix
%
% Output:
% X - estimate of M where we minimize  ||m*ac-Bc||^2 +  ||ar*M-Br||^2
% est_r - estimated rank (when rank r is unknown)
%
function [X, est_r] = SVLS( Br,Bc,Ar,Ac,r, direction )

% Estimate rank using elbow method
if (~exist('r', 'var') || isempty(r))
    fprintf(1,'Rank not specified. Trying to guess ...\n');
    s = svd(Bc);
    [~, r] = max(s(1:end-1)./s(2:end)); est_r = r
end

[uC, ~, ~] = lansvd(Bc,r,'L','OPTIONS'); %find the r largest vectors and singulr values of Bc.
[uR, ~, ~] = lansvd(Br',r,'L','OPTIONS'); %find the r largest vectors and singulr values of Br.

if(~exist('direction', 'var') || isempty(direction))
    direction=0;
end


switch direction
    case 1 % rows
        X=(uR*(Ac'*uR\Bc'))';% find Xr
    case -1 % columns
        X=uC*(Ar*uC\Br); % find Xc
    case 0 % both
        Xr=(uR*(Ac'*uR\Bc'))';% find Xr
        Xc=uC*(Ar*uC\Br); % find Xc
        cost_Xc=norm(Ar*Xc-Br,'fro')+norm(Xc*Ac-Bc,'fro');
        cost_Xr=norm(Ar*Xr-Br,'fro')+norm(Xr*Ac-Bc,'fro');
        
        % Take solution minimizing loss function
        if(cost_Xc<cost_Xr)
            X=Xc;
        else
            X=Xr;
        end
end

