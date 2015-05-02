% Relative Root-Mean-Squared-Error for two matrices
% NEW: Allow computation of squared loss only on subset of indices !! 
% Input:
% X1 - first matrix 
% X2 - second matrix 
% ignore_row_inds - (optional) ignore some rows or entries of matrix when summing
% ignore_col_inds - (optional) ignore some columns of matrix when summing
%
% Output:
% out - RR - || X1-X2||_F / ||X2||_F  (over the set of non-ignored indices)
%
function RR = RRMSE(X1, X2, ignore_row_inds, ignore_col_inds )

if(~exist('ignore_row_inds', 'var') || isempty(ignore_row_inds))
    RR = norm(X1-X2,'fro')/norm(X2,'fro');
else
    if(exist('ignore_col_inds', 'var')) % here input is I and J and we want all X(I,J)
            RR = (norm(X1-X2,'fro') - ...
                sum(sum((X1(ignore_row_inds,ignore_col_inds) - X2(ignore_row_inds,ignore_col_inds)).^2))) / ...
        (norm(X2,'fro') - sum(sum(X2(ignore_row_inds,ignore_col_inds).^2)));
    else % here input is explicitly a list of indices
        if(size(ignore_row_inds, 2)==2) % input is two sets of indices    I(i),J(i)        
            ignore_inds = sub2ind(size(X1), ignore_row_inds(:,1), ignore_row_inds(:,2)); 
        else
            ignore_inds = ignore_row_inds; % input is single 2D indices I(i)
        end
        RR = (norm(X1-X2,'fro') - sum((X1(ignore_inds) - X2(ignore_inds)).^2)) / ...
            (norm(X2,'fro') - sum(X2(ignore_inds).^2));
    end
    
end

