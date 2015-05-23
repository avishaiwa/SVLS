% Sample matrix entries uniformly (matrix-completion model with noise) 
% Input-
% m - random matrix from size nXn
% p-number of measurment that sampled;
% noise- add some normal noise to m_E with std noise.
% Output- 
% m_E - projection of m into the space of matrix with zeros
% everwhere except some random entries.
%
function am = random_measurment(m, p, noise, X_type)
if(~exist('X_type', 'var') || isempty(X_type))
    X_type = 'low_rank';
end
[n1 n2]=size(m); 
am = m + noise*randn(n1,n2);
samp=randsample(n1*n2,n1*n2-p); am(samp) = 0;
switch X_type
    case 'symmetric_low_rank' % symmetrize matrix
        am=am-tril(am)+triu(am)';
end


