%%**********************************************************************
%%  This is a Matlab implementation of an accelerated proximal gradient 
%%  algorithm for solving the regularized semidefinite least squares 
%%  problem:
%%
%%  min_X { f(X) + mu*trace(X) : X psd } 
%%  where f(X) = 0.5*norm(A(X)-b)^2 + 0.5*rho*norm(X,'fro')^2.
%%
%%  At kth iteration, we solve the following subproblem to get Xk:
%%
%%  min_X {0.5*L*norm(X-(Yk-gradf(Yk)/L),'fro')^2 + mu*trace(X) : X psd}
%%  where L is a Lipschitz constant, Yk = beta1*Xk - beta2*Xk_1.
%%
%% Input: nr = size(X,1) 
%%        nc = size(X,2)
%%        Amap = function handle to evaluation A(X)
%%        ATmap = function handle to evaluation the adjoint AT(y)
%%
%% Output: 
%% if X is numeric, then X = standard nrxnc matrix
%% if X is a structure, then matrix = X.U*X.V'; 
%% 
%% NNLS, version 0:  
%% Copyright (c) Feb 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
%%**********************************************************************

   function [X,iter,ttime,sd,runhist] = ...
             APGLsym(nr,nc,Amap,ATmap,bb,mutarget,rhotarget,par)

   clear global  %% important 
   global X Xold Xh Xhold Grad 

   if ~exist('par'); par = []; end
   nmin = min(nr,nc);
   randnstate = randn('state'); randstate = rand('state');
   randn('state',5); rand('state',5);

   tol     = 1e-4;
   maxiter = 100;
   maxrank = min(nmin,100); 
   verbose = 1;
   plotyes = 0;
   linesearch   = 1; 
   continuation = 1; 
   continuation_scaling = 1e-4; 
   truncation     = 1; 
   truncation_gap = 5; 
 
   tstart = clock;
   if isfield(par,'tol'); tol = par.tol; end
   if isfield(par,'maxiter'); maxiter = par.maxiter; end
   if isfield(par,'maxrank'); maxrank = par.maxrank; end
   if isfield(par,'verbose'); verbose = par.verbose; end
   if isfield(par,'plotyes'); plotyes = par.plotyes; end
   if isfield(par,'linesearch');   linesearch = par.linesearch; end   
   if isfield(par,'continuation'); continuation = par.continuation; end
   if isfield(par,'continuation_scaling'); 
      continuation_scaling = par.continuation_scaling; end
   if isfield(par,'truncation');   truncation = par.truncation; end
   if isfield(par,'truncation_gap'); truncation_gap = par.truncation_gap; end
%%
   if ~isa(Amap,'function_handle') | ~isa(ATmap,'function_handle')
      error(['Amap and ATmap must be function handles.']);
   end
%%
%% estimate Lipschitz constant = largest eigenvalue of Amap(ATmap). 
%%
   mm = length(bb); 
   AAT_is_identity = 0; 
   AX = randn(mm,1);
   AY = Amap(ATmap(AX)); 
   if (norm(AX-AY) < 1e-16*norm(AX)); 
      AAT_is_identity = 1; 
      Lipschitz_const = 1; 
   end
   if (AAT_is_identity==0)
      options.tol   = 1e-6; 
      options.issym = true; 
      options.disp  = 0; 
      options.v0    = randn(mm,1);
      Lipschitz_const = eigs(@(y)Amap(ATmap(y)),mm,1,'LM',options); 
      Lipschitz_const = full(Lipschitz_const); 
   end
%%
%% initial mu and target mu
%%
   if (~continuation); continuation_scaling = 1; end
   mu  = mutarget/continuation_scaling;
   rho = rhotarget/continuation_scaling; 
%%
   if (nmin <= 1000)
      fullsvdswitch = 200;
   else
      fullsvdswitch = max(floor(nmin/10)+100,200);
   end
%%
%% initial value of the number of singular values computed
   svp = 1; 
   sv  = svp; 
%%
%% Initialization
%%
   if (max(nr,nc) < 1000);  
      matrix_format = 'standard';
   else
      matrix_format = 'factors'; 
   end   
   if (verbose); 
      fprintf('\n\n  nr = %2.0d, nc = %2.0d, m = %2.0d',nr,nc,mm);
      fprintf('\n  Lipschitz const = %3.2e',Lipschitz_const); 
      fprintf('\n  mu_target = %3.2e',mutarget); 
      fprintf('\n  matrix storage format = %s',matrix_format); 
      fprintf('\n-----------------------------------------------');
      fprintf('-------------------------------------------'); 
      fprintf('\n it     mu       tau      svp    relDist   relRes '); 
      fprintf(' relObjdiff   obj       time')
      fprintf('\n-----------------------------------------------');
      fprintf('-------------------------------------------');  
   end   
   if strcmp(matrix_format,'standard')
      X = sparse(nr,nc);
      Xold = X; 
   else
      Xh.U = sparse(nr,1); Xh.V = sparse(nc,1);
      Xhold = Xh;
   end
   AX = zeros(mm,1); 
   AXold = AX;
   normb = max(norm(bb,'fro'),1); 
   obj   = 0;
   normX = 0;  
   normXold = 0; 
%%
   param.nr = nr; 
   param.nc = nc; 
   param.plotyes = plotyes; 
   param.truncation = truncation;
   param.maxrank = maxrank; 
   param.matrix_format = matrix_format;  
   param.fullsvdswitch = fullsvdswitch; 
%%
   taumax = Lipschitz_const; tau = taumax; taumin = 1e-3*taumax; 
%%
   for iter = 1:maxiter
      if (iter == 1)
         t     = 1;
         beta1 = 1; 
         beta2 = 0; 
         AY = AX; 
      else
         t     = (1+sqrt(1+4*told^2))/2; 
         beta  = (told-1)/t; 
         beta1 = (1+beta)*(1-rho/tau); 
         beta2 = beta*(1-rho/tau); 
         AY    = beta1*AX - beta2*AXold; 
      end
      Gradvec = AY-bb; 
      Grad = ATmap(Gradvec); 
      normGrad = norm(Grad,'fro'); 
      %%
      %% G = Y-Grad/tau, where Y = beta1*X-beta2*Xold. 
      %%
      if strcmp(matrix_format,'standard')
         trXXold = sum(sum(X.*Xold));     
      else
         trXXold = sum(sum((Xhold.U'*Xh.U).*(Xhold.V'*Xh.V)));
      end
      trXGrad    = AX'*Gradvec;
      trXoldGrad = AXold'*Gradvec;
      normY2  = beta1^2*normX^2 + beta2^2*normXold^2 - 2*beta1*beta2*trXXold;
      trYGrad = beta1*trXGrad - beta2*trXoldGrad;
      normG2  = normY2 - (2/tau)*trYGrad + (1/tau^2)*normGrad^2;  
      %%  
      %% preditermined number of singular values
      %%
      param.beta1 = beta1; param.beta2 = beta2; param.beta3 = 1/tau;
      if (iter>1); param.verbose=verbose; else; param.verbose=0; end
      if (svp == sv) & (mu/mutarget < 50)
         sv = min(svp+5,nmin); 
      else
         sv = min(svp+1,nmin);
      end
      if (mu > mutarget);
         svdtol = min(1e-3,0.01*mu);  
      else
         svdtol = min(0.1*tol,0.01*mu); 
      end
      if (iter <= 10); 
         gap = inf; 
      else
         gap = truncation_gap; 
      end
      [U,V,sd,svp] = proxmapsym(sv,1,mu/tau,svdtol,gap,param);
      %%
      if strcmp(matrix_format,'standard')
         trXnewX    = sum(sum(U.*(X*V))); 
         trXnewXold = sum(sum(U.*(Xold*V)));
      else
         trXnewX    = sum(sum((Xh.U'*U).*(Xh.V'*V)));
         trXnewXold = sum(sum((Xhold.U'*U).*(Xhold.V'*V)));
      end
      told     = t;      
      normXold = normX; 
      objold   = obj;    
      AXold    = AX;
      %% Note: must be after AXold = AX; 
      if strcmp(matrix_format,'standard')
         Xnew  = U*V';
         AX = Amap(Xnew); 
      else
         Xhnew.U = U; Xhnew.V = V; 
         AX = Amap(Xhnew); 
      end
      trXnewGrad = AX'*Gradvec; 
      trXnewY = beta1*trXnewX - beta2*trXnewXold; 
      trXnewG = trXnewY - (1/tau)*trXnewGrad;      
      normXG2 = norm(sd)^2 + normG2 - 2*trXnewG; 
      objlinesearch = 0.5*normGrad^2 - normGrad^2/(2*tau) + (tau/2)*normXG2;  
      objlinesearch = objlinesearch + 0.5*rho*normY2 + mu*sum(sd); 
      %%
      if strcmp(matrix_format,'standard')
         Xold = X; X = Xnew;
      else
         Xhold = Xh; Xh = Xhnew; 
      end
      trXXold   = trXnewX; %% important
      normX     = norm(sd,'fro');
      normRes   = norm(AX-bb,'fro'); 
      obj       = 0.5*normRes^2 + mu*sum(sd) + 0.5*rho*normX^2;
      normXdiff = sqrt(normX^2 + normXold^2 - 2*trXXold);
      normX_Y   = sqrt(normX^2 +normY2 -2*trXnewY);
      if (AAT_is_identity)
         normS = (1/tau)*sqrt(tau^2*normX_Y^2 + (1-2*tau)*norm(AX-AY)^2);
      else
         AX_AY = AX-AY; 
         normtmp = norm(ATmap(AX_AY),'fro'); 
         normS = (1/tau)*sqrt(tau^2*normX_Y^2 -2*tau*norm(AX_AY)^2 +normtmp^2);
      end
      relDist    = normS/max(normX,1);
      relObjdiff = abs(obj-objold)/max(1,obj);  
      relRes     = normRes/normb; 
      runhist.obj(iter)           = obj; 
      runhist.objlinesearch(iter) = objlinesearch; 
      runhist.relDist(iter)      = relDist;
      runhist.relObjdiff(iter)    = relObjdiff;  
      runhist.relRes(iter)        = relRes;  
      runhist.tau(iter)           = tau;
      runhist.mu(iter)            = mu; 
      runhist.svp(iter)           = svp; 
      runhist.sv(iter)            = sv; 
      %%
      if (verbose)
         fprintf('\n %2.0d  %3.2e %3.2e |%3.0d %3.0d| ',iter,mu,tau,svp,sv);
	 fprintf('%3.2e %3.2e %3.2e| %5.4e|',relDist,relRes,relObjdiff,obj); 
      end    	  
      %% 
      %% check stopping criterion
      %%
      if (mu == mutarget) & (iter > 4)
         if (max([relDist, 0.1*relObjdiff]) < tol) 
	    msg = sprintf('[a] relDist < %3.2e',tol);
 	    fprintf('\n %s',msg); 
            break; 
         end
	 if (max([0.5*relDist, 100*relObjdiff]) < tol)  
	    msg = sprintf('[b] relObjdiff < %3.2e',0.01*tol);
 	    fprintf('\n %s',msg); 
            break; 
         end
         idx = [iter-2:iter];
         relResratio = runhist.relRes(idx-1)./runhist.relRes(idx); 
         if (all(abs(relResratio-1) < 5*tol) & (relObjdiff < 0.1*tol)) ...
            & (relDist < 10*tol)
	    msg = sprintf('[c] lack of progress in relRes');
            fprintf('\n %s : %2.1e',msg,max(relResratio-1)); 
            break; 
         end
      end
      %%
      %% perform linesearch to update tau. 
      %% 
      if (mu < 100*mutarget) & (linesearch)
         if (verbose); fprintf(' %5.4e|',objlinesearch); end
         const_eta = 0.8; 
         if (obj < objlinesearch) 
   	    tau = min(taumax,max(const_eta*tau,taumin)); 
         else
            %% restart using old X. 
            if strcmp(matrix_format,'standard')
               X = Xold;   
            else
               Xh = Xhold;   
            end
            AX     = AXold; 
            normX  = normXold;
            if (mu == mutarget); 
               taumin = tau/const_eta; 
            end; 
            told = 1; 
            tau  = min(tau/const_eta,taumax); 
         end
         %% not a good idea, can cause strange error for RKE problem
         %%if (mu == mutarget) & (iter > 4)
         %%   idx = [iter-4:iter]; 
         %%   if all(runhist.obj(idx) < 0.98*runhist.objlinesearch(idx))
         %%      taumin = taumin*const_eta; 
         %%   end
         %%end
      end
      %%    
      %% update mu
      %%
      mu = max(0.7*mu,mutarget); 
      if (verbose); fprintf(' %3.2e|',etime(clock,tstart)); end
   end 
   sinmax = max(sd);
   sinmin = min(sd(find(sd>0)));
   ttime  = etime(clock,tstart);
   if ~strcmp(matrix_format,'standard')
      X = Xh; 
   end
   randn('state',randnstate); rand('state',randstate);
%%
%% Print results
%%
   if (verbose)
      fprintf(1,'\n Finished the main algorithm!\n')
      fprintf(1,' Objective function        = %6.5e\n',obj);
      fprintf(1,' Number of iterations      = %2.0d\n',iter);
      fprintf(1,' Number of singular values = %2.0d\n',svp);
      fprintf(1,' max singular value        = %3.2e\n',sinmax);
      fprintf(1,' min singular value        = %3.2e\n',sinmin);
      fprintf(1,' CPU time                  = %3.2e\n',ttime);
      fprintf(1,' norm(X-Xold,''fro'')/norm(X,''fro'') = %3.2e\n',relDist);
      fprintf(1,'\n');
   end
%%**********************************************************************

