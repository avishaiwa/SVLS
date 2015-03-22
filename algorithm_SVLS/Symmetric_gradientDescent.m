% Minimize f(x,y) = ||Ar*X*Y'-Br||^2 using the
% gradient-descent algorithm
%
% Input:
% X_hat - starting point of the algorithm (product U*U' for starting point)
% B = A*X + Normal(0,noise)
% A = matrix size kXn - row measurements
% r -  rank r of unknown matrix N
% alpha - gradient step size
% tol - tolerance for stopping criteria
% max_iter - maximum # of iterations allowed
%
% Output
% N - estimation of m without any noise
% dis - matrix of distances as function of iteration
%
function [X,dis] = Symmetric_gradientDescent( X_hat, B, A, r, tol, max_iter )
% using steepest decsent to estimate X when X*ac = Bc and ar*X=Br (whithout noise)
addpath('PROPACKmod')

if(~exist('tol', 'var') || isempty(tol))
    tol = 0.00001; % set default tolerance
end
if(~exist('max_iter', 'var') || isempty(max_iter))
    max_iter = 100; % set default # of iterations
end


[u, d, v] = lansvd(X_hat,r,'L','OPTIONS'); % find the r largest vectors and singular values.
u = u*sqrt(d);
v = v*sqrt(d);
optT = 0.5*(u+v); % Temp. We want to find the closest positive-definite rank-r approximation for N_hat. Is there a closed-form solution for this?
% (current heuristic: perform S.V.D. and take average of u and v)
norm(optT)

X1 = X_hat;
dis=zeros(1,max_iter);
iter=0;
X0 = zeros(size(X1));

% gradient descent
while(norm(X1-X0,'fro')>tol&&iter<max_iter)
    w = Symmetric_gradLs(optT,B,A);
    norm(X1-X0,'fro')
    alpha = getoptAlpha(T,W,S,Z,B,A);
    optT = optT + alpha*w; %gradient step
    X0 = X1;
    X1 = optT*optT';
    iter = iter+1;
    dis(iter)=norm(A*X1-B,'fro');
end
X = optT*optT';

end

% Internal function: the gradient of the cost function.
function W = Symmetric_gradLs(T,B,A)
%W = ar'*(ar*T*S-Br)*S'+(T*S*ac-Bc)*ac'*S';

% Compute symmetric gradient:
% Grad(L(U ; A, B) ) =  2 [ A^T A U U^T + U U^T A^T A - A^T B - B^T A ] U
W = (A'*A*(T*T') + T*T'*(A'*A) - A'*B - B'*A )* T;


end

%loss function
function out = F_t(T,S,B,A)

out=0.5*norm(A*T*S-B,'fro')^2;
end

%step size
function out = getoptAlpha(T,W,S,Z,B,A)
norm2WZ = norm(W,'fro')^2 + norm(Z,'fro')^2;
f(1) = F_t(T,S,B,A);

t = -1e-1 ;
for i = 1:20
    f(i+1) = F_t(T+t*W,S+t*Z,B,A) ;
    
    if( f(i+1) - f(1) <= .5*(t)*norm2WZ )
        out = t ;
        return;
    end
    t = t/2 ;
end
out = t ;
end




