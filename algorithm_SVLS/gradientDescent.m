% Minimize f(x,y) = ||Ar*W*S'-Br||^2 + ||W*S'*Ac-Bc|^2|  using the
% gradient-descent algorithm.
% Note: this assumes that elements which were sample twice (in their row and column)
% are considered twice in the loss function.
%
% Input:
% X_hat - starting point of the algorithm (product X*Y' for starting point)
% Br = Ar*X + Normal(0,noise)
% Bc = X*Ac + Normal(0,noise)
% Ar = matrix size pXn - row measurements
% Ac = matrix size nXp - column measurements
% r -  rank r of unknown matrix N
% alg_str - which algorithm to use (default is standard gradient descent)
% tol - tolerance for stopping criteria
% max_iter - maximum # of iterations allowed
%
% Output:
% X - estimation of m without any noise
% dis - matrix of distances (errors, divided into two terms) as function of iteration
%
function [X,dis] = gradientDescent( X_hat,Br,Bc,Ar,Ac,r, ...
    alg_str,tol, max_iter )
% using steepest decsent to estimate X when X*ac = Bc and ar*X=Br (whithout noise)
if(~exist('tol', 'var') || isempty(tol))
    tol = 0.00001; % set default tolerance
end
if(~exist('max_iter', 'var') || isempty(max_iter))
    max_iter = 100; % set default # of iterations
end
if(~exist('alg_str', 'var') || isempty(alg_str))
    alg_str = 'gradient_descent'; % set default algorithm
end

% Problem! for small n lansvd sometimes fails! replace it with 'standard' SVD
if(size(X_hat,1) < 100)
    [u, d, v ] = svd(X_hat); u=u(:,1:r); v=v(:,1:r); d=d(1:r,1:r);
else
    [u, d, v] = lansvd(X_hat,r,'L','OPTIONS'); %find the r largest vectors and singular values.
end
u = u*sqrt(d);
v = v*sqrt(d);
optT = u;
optS = v';
if(norm(u*v' - X_hat, 'fro')/norm(X_hat) > 0.1)
    error_in_propack = 99999
end

X1 = X_hat;
dis=zeros(1,max_iter);
iter=0;
X0 = randn(size(X1));
k_r = size(Ar,1); k_c = size(Ac,2); % get number of measurements

% Gradient descent to improve solution
while( (norm(X1-X0,'fro')>tol) && (iter<max_iter))
    switch alg_str
        case 'gradient_descent'               
            [w, z] = GradLs(optT,optS,Br,Bc,Ar,Ac);
        case 'stochastic_gradient_descent'
            [w, z] = StocGradLs(optT,optS,Br,Bc,Ar,Ac,randi(k_r),randi(k_c));
    end
    alpha = getoptAlpha(optT,w,optS,z,Br,Bc,Ar,Ac); %%line search for the best gradient step
    optT = optT + alpha*w; %gradient step
    optS = optS + alpha*z; %gradient step
    X0 =X1;
    X1 = optT*optS;
    iter = iter+1;
    fprintf('norm(X-X_old)=')
       norm(X1-X0,'fro')

    %%loss function
    dis(iter)=F_t(optT,optS,Br,Bc,Ar,Ac);
end
X = optT*optS;

end

% Internal function: the gradient of the cost function.
function [W, Z] = GradLs(T,S,Br,Bc,ar,ac)
W = ar'*(ar*T*S-Br)*S'+(T*S*ac-Bc)*ac'*S';
Z = T'*ar'*(ar*T*S-Br)+T'*(T*S*ac-Bc)*ac';
end

% Internal function: stochastic gradient of the cost function,
% considering only one row (i-th) and one column (j-th)
function [W, Z] = StocGradLs(T,S,Br,Bc,ar,ac,i,j)
W = ar(i,:)'*((ar(i,:)*T)*S-Br(i,:))*S'+(T*(S*ac(:,j))-Bc(:,j))*ac(:,j)'*S';
Z = T'*ar(i,:)'*((ar(i,:)*T)*S-Br(i,:))+T'*(T*(S*ac(:,j))-Bc(:,j))*ac(:,j)';
end

%loss function
function out = F_t(T,S,Br,Bc,Ar,Ac)

out=0.5*norm(Ar*T*S-Br,'fro')^2+0.5*norm(T*S*Ac-Bc,'fro')^2;
end

%Compute optimal step size - based on ???
function out = getoptAlpha(T,W,S,Z,Br,Bc,Ar,Ac)
norm2WZ = norm(W,'fro')^2 + norm(Z,'fro')^2;
f(1) = F_t(T,S,Br,Bc,Ar,Ac);

t = -1e-1 ;
for i = 1:20
    f(i+1) = F_t(T+t*W,S+t*Z,Br,Bc,Ar,Ac) ;
    
    if( f(i+1) - f(1) <= .5*(t)*norm2WZ )
        out = t ;
        return;
    end
    t = t/2 ;
end
out = t ;
end



