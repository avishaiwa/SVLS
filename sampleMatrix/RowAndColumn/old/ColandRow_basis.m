% Input:
% n - number of measurements
% r - rank of unknown matrix
% p - number of variables
% noise - standard deviation of Gaussian additive noise
%
% Output:
%
% m - matrix n*n from rank r
% Br = ar*m + Normal(0,noise)
% Bc = m*ac + Normal(0,noise)
% ar = matrix size pXn - single entries row measurements
% ac = matrix size nXp - single entries column measurements
%
function [ m,Br,Bc,ar,ac ] =ColandRow_basis( n,r,p,noise, X_type)

n1=0;n2=0;
if length(n)==2
    n1 =n(1);n2=n(2);
elseif  length(n)==1
    n1=n;n2=n;
end
        

if(~exist('X_type', 'var') || isempty(X_type))
    X_type = 'low_rank';
end

u = randn(n1,r)/r^0.25;
v = randn(n2,r)/r^0.25;

switch X_type
    case 'low_rank'
       m = u*v';
    case 'symmetric_low_rank'
       m = u*u'; 
end

am = m +noise*randn(n1,n2);
row_samp=randsample(n1,p);
ar = zeros(p,n1);
ac = zeros(n2,p);
for i=1:p
        ar(i,row_samp(i))=1;
end
col_samp=randsample(n,p);
for i=1:p
        ac(col_samp(i),i)=1;
end
Br = ar*am;
Bc =am*ac;

switch X_type
    case 'low_rank'
       m = u*v';
    case 'symmetric_low_rank'
       ar=ac';
       Br=Bc';
end



end


