% Minimize F(W,S) = ||Ar*W*S'-Br||_F^2 + ||W*S'*Ac-Bc||_F^2|  using the
% alternating least-squares algorithm.
% Note: this assumes that elements which were sampled twice (in their row and column)
% are considered twice in the loss function.
%
% Input:
% X_hat - starting point of the algorithm (can use product X*Y' for starting point)
% Br = Ar*X + Normal(0,noise)
% Bc = X*Ac + Normal(0,noise)
% Ar = matrix size pXn - row measurements
% Ac = matrix size nXp - column measurements
% r -  rank r of unknown matrix N
% tol - tolerance for stopping criteria
% max_iter - maximum # of iterations allowed
%
% Output:
% X - estimation of m without any noise
% dis - matrix of distances (errors, divided into two terms) as function of iteration
%
function [X, dis] = Alternating_LS( X_hat, Br, Bc, Ar, Ac, r, ...
    tol, max_iter)

% using steepest decsent to estimate X when X*ac = Bc and ar*X=Br (whithout noise)
if(~exist('tol', 'var') || isempty(tol))
    tol = 0.00001; % set default tolerance
end
if(~exist('max_iter', 'var') || isempty(max_iter))
    max_iter = 100; % set default # of iterations
end

% Problem! for small n lansvd sometimes fails! replace it with 'standard' SVD
if(size(X_hat,1) < 100)
    [U, d, V ] = svd(X_hat); U=U(:,1:r); V=V(:,1:r); d=d(1:r,1:r);
else
    [U, d, V] = lansvd(X_hat,r,'L','OPTIONS'); %find the r largest vectors and singular values.
end
U = U*sqrt(d);
V = V*sqrt(d);

% Check for decomposition
if(norm(U*V' - X_hat, 'fro')/norm(X_hat) > 0.1)
    error_in_propack = 99999
end

prev_X = X_hat;
dis=zeros(1,max_iter);
iter=0;
new_X = randn(size(prev_X)); % why take a random matrix?
[k_r, n1] = size(Ar); [n2, k_c] = size(Ac); % get number of measurements


% Create big matrices with Kronecker products
% big_A = kron(Ar, eye(k_r)) + kron(Ac', eye(k_c)); % make big measurements matrix
% big_b = [mat2vec(Br) mat2vec(Bc')]; % vectorize measurements
[big_A, big_b] = row_column_to_general_affine(Ar, Ac, Br, Bc);

% perform alternating minimization to improve solution
while( (norm(prev_X-new_X,'fro')>tol) && (iter<max_iter))
    [U_mat, V_mat, x_vec2] = bilinear_to_general_affine(U, V); % compute representation
    
    new_V = (big_A*U_mat) \ big_b'; old_V=V; V=vec2mat(new_V, r); % change V while U is kept fixed
    
    run_kronecker_least_squares=1; 
    if(run_kronecker_least_squares) % NEW! try to simplfy U matrix using Kronecker products 
%        big_A_times_U_mat2 = [kron(eye(n2), (Ar*U)') kron(Ac, U')]';
        big_A_times_U_mat2 = [kron(eye(n2), (Ar*U)) ; kron(Ac', U)];
        A_trans_A = big_A_times_U_mat2'*big_A_times_U_mat2; 
        A_trans_A2 = kron(eye(n2), (Ar*U)'*(Ar*U)) + kron(Ac*Ac', U'*U);            
        A_trans_A_inv = inv(A_trans_A);
    
        % Compute generalized SVD (take transpose) 
        [U11,U12,V1,Diag11,Diag12] = gsvd(eye(n2),Ac');
        [U21,U22,V2,Diag21,Diag22] = gsvd((Ar*U),U);
        
        should_be_zero = max(max(abs(U11*Diag11*V1' - eye(n2) ) ) )
        should_be_zero = max(max(abs(U12*Diag12*V1' - Ac' ) ) )
        should_be_zero = max(max(abs(U21*Diag21*V2' - Ar*U ) ) )
        should_be_zero = max(max(abs(U22*Diag22*V2' - U ) ) )
        
        % Next, perform QR-decomposition for the matrices: 
        % Compute three matrices composing A 
        % Wrong: from paper: 
        % A_part1 = [kron(U11, U12) zeros(size(U11,1)*size(U12,1), size(U21, 2)*size(U22, 2)) ; ...
        %    zeros(size(U21,1)*size(U22,1), size(U11, 2)*size(U12, 2)) kron(U21, U22)]; 
        % A_part2 = [kron(Diag11, Diag12)  kron(Diag21, Diag22)]'; 
        % My version: 
        A_part1 = [kron(U11, U21) zeros(size(U11,1)*size(U21,1), size(U12, 2)*size(U22, 2)) ; ...
            zeros(size(U12,1)*size(U22,1), size(U11, 2)*size(U21, 2)) kron(U12, U22)]; 
        A_part2 = [kron(Diag11, Diag21) ; kron(Diag12, Diag22)]; 
        
        A_part3 = kron(V1', V2'); 
        big_A_times_U_mat3 = A_part1*A_part2*A_part3; 
        should_be_zero_big = max(max(abs(big_A_times_U_mat2 - big_A_times_U_mat3 ) ) )

        
    end % if run_kronecker_least_squares
    
    intermediate_dis = row_column_l2_loss(U,V',Br,Bc,Ar,Ac)
    [U_mat, V_mat, x_vec2] = bilinear_to_general_affine(U, V); % compute representation again
    new_U = (big_A*V_mat) \ big_b';  old_U=U; U=vec2mat(new_U, size(X_hat,1))'; % change U while V is kept fixed
    
    prev_X = new_X;
    new_X = U*V';
    fprintf('norm(X-X_old)=%3.3f\n', norm(new_X-prev_X,'fro')) % why does this keep changing? 
    iter = iter+1;
    
    % Compute loss function
    dis(iter)=row_column_l2_loss(U,V',Br,Bc,Ar,Ac);
    fprintf('Iter=%d L2 Loss=%3.3f\n', iter, dis(iter))
end
X = U*V'; % compute last time

dis=dis(1:iter); % get number of iterations 

end % function




