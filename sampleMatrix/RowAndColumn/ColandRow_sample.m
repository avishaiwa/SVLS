%create matrix with column and row samples for optSpace
% Input:
% m - true matrix from size nXn
% p - number of rows and columns that sampeled
% noise- add some normal noise to m_E with std noise.
% X_type - type of matrix m
%
% Output: 
% m_E - m_E is the projection of m into the space of matrix with zeros
% everwherhe accept some randoms rows and columns.
%
function m_E = ColandRow_sample(m, k, noise, X_type)

if(~exist('X_type', 'var') || isempty(X_type))
    X_type = 'low_rank';
end
[n1 n2]=size(m); 
row_samp=sort(randsample(n1,k));
switch X_type
    case 'low_rank'
        col_samp=sort(randsample(n2,k));
    case 'symmetric_low_rank'
        col_samp=row_samp; 
end

% Sample entries at row and columns
m_E = zeros(n1,n2);
m_E(row_samp,:)=m(row_samp,:) + noise.*randn(k,n2);
m_E(:,col_samp)=m(:,col_samp) + noise.*randn(n1,k);

