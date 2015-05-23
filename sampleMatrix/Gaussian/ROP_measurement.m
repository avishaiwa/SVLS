% Sample Rank-One-Projection noisy measurements of a matrix  
% beta^T X gamma for Gaussian i.i.d. Gaussian vectors beta, gamma
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
% b = vector of measurements: b(i) = Ar(i,:) X Ac(:,i) + Normal(0,noise)
% Ar = matrix size pXn - row vectors for projections
% Ac = matrix size nXp - column vectors for projections
%
function [b, Ar, Ac] = ROP_measurement(X, k, noise, X_type)

% Set default parameters
if(~exist('X_type', 'var') || isempty(X_type))
    X_type = 'low_rank';
end
n=length(X);

% Choose measurements
Ac = randn(n,k); b=zeros(k,1); 
switch X_type
    case 'low_rank'
        Ar = randn(k,n);
    case  'symmetric_low_rank'
        Ar = Ac';  % force Ar=Ac' in symmetric case
end
for i=1:k % calculate projections (can't avoid a loop here?)  
    b(i) = Ar(i,:)*X*Ac(:,i); 
end
b=b+noise*randn(k,1); 


