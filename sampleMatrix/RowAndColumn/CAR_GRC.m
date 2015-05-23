% Simulation of matrix noisy measurements (columns and rows with Gaussian coefficients)
%
% Input:
% X - 'true' square matrix n*n
% k - number of measurements
% noise - standard deviation of Gaussian additive noise
% X_type - (NEW) type of matrix to simulate. Default is just an mxn matrix.
%           Can allow 'symmetric'
% measurement_type - (NEW) type of measurements to perform
%
% Output:
%
% Br = ar*X + Normal(0,noise)
% Bc = X*ac + Normal(0,noise)
% Ar = matrix size pXn - row measurements
% Ac = matrix size nXp - column measurements
%
function [Br, Bc, Ac, Ar] = CAR_GRC(X, k, noise, X_type)

% Set default parameters
if(~exist('X_type', 'var') || isempty(X_type))
    X_type = 'low_rank';
end
[n1 n2]=size(X);

% Choose measurements
Ac = randn(n2,k)/(n2^0.5) ;%/(n1^0.25) ;
Bc = X*Ac+ noise*randn(n1,k);
switch X_type
    case 'low_rank'
        Ar = randn(k,n1)/(n2^0.5) ;%/(n1^0.25) ;
        Br = Ar*X+ noise*randn(k,n2); % Create noisy measurements
    case  'symmetric_low_rank'
        Ar = Ac';  % force Ar=Ac' in symmetric case
        Br=Bc';
end




