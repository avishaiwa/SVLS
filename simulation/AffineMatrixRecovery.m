% Perform one simulation
%
% Input:
% X - the true matrix (optional, not used by algorithm - just for error evaluation)
% am - measurements
% alg - String - the algorithm to implement
% r - the rank of m
% k - number of row and columns that are sampled
% noise - STD of the measurements
% maximum_iter - maximum iteretion in gradient descent
% tol - stopping rule for greadient descent
% measure - affine measurements
% r - rank (optional) ??
%
% Output:
% X_hat - estimated matrix
% time - number of second the algorithm alg ran
% rel_error - RRMSE(X, hat_X) (if X is given in input)
% iter_num - number of iterations
%
function [X_hat, time, rel_error, iter_num] = ...
    AffineMatrixRecovery(X, alg, k, noise, maximum_iter, tol, measure)
[n1, n2]=size(X);
r = rank(X);
tic;
switch lower(alg) % Choose algorithm for reconstruction
    
    %%SVLS algorithm         %SVLS and gradient descent with SVLS as starting point, stochastic gradient descent
    case {'svls', 'gradient_svls', 'stochastic_gradient_svls'} % Our algorithm
        [X_hat, r] = SVLS( measure.Br, measure.Bc, measure.Ar, measure.Ac,[]);
        iter_num=1 ;
        
        %%iterative SVLS algorithm
    case 'iterative_svls' % Our algorithm
        [X_hat, ~, row_loss_vec, col_loss_vec] = ...
            Iterative_SVLS(measure.Br, measure.Bc, measure.Ar, measure.Ac, [], 10); % what is r???
        
        %figure; hold on;
        %plot(row_loss_vec); plot(col_loss_vec, 'g'); plot(row_loss_vec+col_loss_vec, 'r');
        %legend({'row-loss', 'col-loss', 'Combined'}); xlabel('Iteration'); ylabel('Loss');
        
        iter_num=1 ;
        
        
    case 'optspacenotrim'         %%optSpace without the trimming step
        [W, S, T, iter_num] = OptSpacenoTrim(measure.am,r,maximum_iter,tol); %normal optSpace
        X_hat= W*S*T';
        iter_num = length(iter_num);
        
        %%minimize the loss function
        %%(||Ar(X)-Br||_Fro)^2+(|(X)Ac-Bc||_Fro)^2 with zeros as starting point
    case 'gradient_random_start'
        X_hat =randn(n1,n2);
       
        %%SVLS and gradient descent on symmetric matrix
    case  'symmetric_gradient_descent'
        [N, r] = SVLS( measure.Br, measure.Bc, measure.Ar, measure.Ac,[]);
        B=0.5*(N+N'); % symmetrize matrix
        N = 0.5*(B+real(sqrtm(B'*B)));
        [X_hat, iter_num] = Symmetric_gradientDescent( N, measure.Br, measure.Ar, r, tol, maximum_iter );%gradient descent with the putput of basis as starting point
        iter_num= length(iter_num);
        
        %%SVLS and stochastic gradient descent
       
        
        %%SVT algorithm
    case 'svt' % Singular Value thresholding
        % Set default parameters, The parmater is taken as in the papaer A Singular Value
        %Thresholding Algorithm for Matrix Completion
        tau = 5*sqrt(n1*n2);
        k1 = (n1+n2)*k-k^2;
        delta = 1.2*n1*n2/k1;
        [X_tmp,S,Y_tmp,numiter,out]  = SVT([n1 n2], find(measure.am), measure.am(measure.am~=0), tau,delta,maximum_iter,tol,10^-7); %svt algorithm
        X_hat= X_tmp*S*Y_tmp';
        iter_num = numiter;
        
         case 'svt_col_row' % Singular Value thresholding
        % Set default parameters, The parmater is taken as in the papaer A Singular Value
        %Thresholding Algorithm for Matrix Completion
        tau = 5*sqrt(n1*n2);
        k1 = (n1+n2)*k-k^2;
        delta = 1;
        [X_tmp,S,Y_tmp,numiter,out]  = SVT([n1 n2], find(measure.am), measure.am(measure.am~=0), tau,delta,maximum_iter,tol,noise+10^-7); %svt algorithm
        X_hat= X_tmp*S*Y_tmp';
        iter_num = numiter;
        
        %%optSpace algorithm
    case 'optspace' % The optSpace algorithm from ..
        [X_tmp, S, Y_tmp, iter_num] = OptSpace(measure.am,[],maximum_iter,tol); %normal optSpace
        X_hat= X_tmp*S*Y_tmp';
        iter_num = length(iter_num);
        
        
        %%APGL algorithm
    case 'apgl' % The APGL algorithm from ..
        par.maxiter=maximum_iter;
        par.tol = tol;
        tau=noise*10;
        [X_hat,numiter,~,~,~] =  APGL(n1,n2,'NNLS', measure.map, measure.tmap, measure.bb,tau,0,par); % minimize nuclear norm
        iter_num= numiter;
        
end % switch algorithm
RRMSE(X_hat,X)
% apply a post-processing step
if(strfind(lower(alg), 'gradient')) % Choose algorithm for reconstruction % Next step: perform gradient descent to improve solution
    if(~exist('X_hat', 'var'))
        X_hat = X;
    end
    
     grad_alg = str2word('_', lower(alg), 1)%%First word in the algorithm

     %If the first word is gradient make gradient descent and if the first
     %woed is stotchastic implement SGD
    if strcmp(grad_alg,'gradient')
        grad_alg =  'gradient_descent';
    end
     if strcmp(grad_alg, 'stochastic')
        grad_alg =  'stochastic_gradient_descent';
     end
    
    [X_hat, iter_num] = gradientDescent(X_hat, measure.Br, measure.Bc, measure.Ar, measure.Ac, r, ...
        grad_alg, tol, maximum_iter); %gradient descent with the putput of basis as starting point
    iter_num= length(iter_num);
end

time = toc;
%score = norm( measure.Ar*X_hat-measure.Br,'fro')+norm( X_hatmeasure.Ar-measure.Bc,'fro')
if(exist('X', 'var') && (~isempty(X)))
    rel_error =  RRMSE(X_hat,X); % compute estimation error
else
    rel_error = [];
end

