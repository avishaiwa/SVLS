%%************************************************************************
%% run multivariate linear regression problems: 
%% 
%% min { 0.5*norm(H-G*X,'fro')^2 + mu*sum(svd) } 
%% where G (nrxnr), H (nrxnc) are given matrices
%%
%% For details, see: 
%% [1] M. Yuan, A. Ekici, Z. Lu and R.D.C. Monteiro,  
%%     Dimension Reduction and Coefficient Estimation in the 
%%     Multivariate Linear Regression, 
%%     Journal of the Royal Statistical Society, Series B, 69(3):329-346, 2007.
%%
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
%%************************************************************************

    clear all;

    eg = 2; 
    if (eg==1) 
       nr = 200; nc = 200; 
       p = nr;  
       G = diag(rand(nr,1)); 
       rr = 20;  
    elseif (eg==2)
       nr = 500; nc = 250; 
       p = 2500; 
       G = rand(p,nr); G = G/norm(G,'fro');
       rr = 50; 
    end
    M.U = randn(nr,rr); M.V = randn(nc,rr); 
    normM = sqrt(sum(sum((M.U'*M.U).*(M.V'*M.V))));
    randmat = randn(p,nc);
    H = (G*M.U)*M.V'; 
    H = H + 0.1*norm(H,'fro')*randmat/norm(randmat,'fro');
    isdiagG = (p==nr) & (norm(G-spdiags(diag(G),0,p,nr),'fro')==0); 
    if ~isdiagG
       [U,S,V] = mexsvd(G,0); 
       G = S; H = U'*H; 
       M.U = V'*M.U; 
    end
%%
    tstart = clock;
    bb = H(:); 
    Amap  = @(X) Amap_MLR(X,G);
    ATmap = @(y) ATmap_MLR(y,G);
    B = ATmap(bb); 
    options.tol = 1e-8; 
    mumax = svds(sparse(B),1,'L',options);
    mu_scaling = 1e-4;  
    mutarget   = mu_scaling*mumax;
    par.tol     = 1e-4;
    par.verbose = 1;
    par.plotyes = 1;
    par.continuation_scaling = mu_scaling;  
    par.truncation       = 1; 
    par.truncation_gap   = 20; 
    par.maxiter  = 200; 
    problem_type = 'NNLS'; 
    [X,iter,time,sd,hist] = ...
     APGL(nr,nc,problem_type,Amap,ATmap,bb,mutarget,0,par);
%%
    if exist('M'); 
       if isstruct(X)
          normX = sqrt(sum(sum((X.U'*X.U).*(X.V'*X.V))));
          trXM = sum(sum((M.U'*X.U).*(M.V'*X.V)));
       else
          normX = norm(X,'fro'); trXM = sum(sum(M.U.*(X*M.V))); 
       end
       mse = sqrt(normX^2+normM^2-2*trXM)/normM;
       runhist.mse   = mse; 
    end
    runhist.time  = etime(clock,tstart);
    runhist.mu    = mutarget; 
    runhist.mumax = mumax; 
    runhist.iter  = iter;
    runhist.obj   = hist.obj(end);
    runhist.svp   = hist.svp(end);
    %%
    %% report results in a table
    %%
    fprintf(' Problem: nr = %d, nc = %d',nr,nc); 
    fprintf('\n-----------------------------------------------');
    fprintf('------------------------------')
    fprintf('\n problem type       :  %s',problem_type);
    fprintf('\n iterations         :  %6.0f',runhist.iter);
    fprintf('\n # singular         :  %6.0f',runhist.svp);
    fprintf('\n obj  value         :  %6.5e',runhist.obj);
    fprintf('\n cpu   time         :  %6.2e',runhist.time);
    fprintf('\n mean square error  :  %6.2e',runhist.mse);
    fprintf('\n------------------------------------------------'); 
    fprintf('------------------------------\n')
%%************************************************************************
