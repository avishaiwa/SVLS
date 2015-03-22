%Sample Gaussian affine measurements of m map(m) = am  when map is a matrix a
%and evey entry in a is N(0,1).
% Input:
% m - matrix of size nXn
%
% Output
% map - matrix from size pXn^2 as described above map(m) = a*vec(m).
% tmap = a'  (transpose of map)
% am = map(m) + noise*N(0,1)
%
function [map, tmap, am] = Normal_measurment(m,n,p,noise)
n1=0;n2=0;
if length(n)==2
    n1 =n(1);n2=n(2);
elseif  length(n)==1
    n1=n;n2=n;
end
        
    
a = randn(p,n1*n2)/sqrt(n1*n2); % draw coefficients 
am = a*m(:)+noise*randn(p,1);
map = @(x)a*x(:); % represent a as a transformation (not matrix) 
tmap = @(x)reshape(a'*x,n1,n2);

