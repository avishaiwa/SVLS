%%**********************************************************************
%%  This is a Matlab implementation of an accelerated proximal gradient 
%%  algorithm for solving the nuclear norm regularized least squares 
%%  or regularized semidefinite least squares problem: 
%%
%%  min_X {f(X) + mu*sum(svd(X))}
%%  or 
%%  min_X { f(X) + mu*trace(X) : X psd } 
%%  where f(X) = 0.5*norm(A(X)-b)^2 + 0.5*rho*norm(X,'fro')^2.
%%
%%  At kth iteration, we solve the following subproblem to get Xk:
%%
%%  min_X {0.5*L*norm(X-(Yk-gradf(Yk)/L),'fro')^2 + mu*sum(svd(X))}
%%  or
%%  min_X {0.5*L*norm(X-(Yk-gradf(Yk)/L),'fro')^2 + mu*trace(X) : X psd}
%%  where L is a Lipschitz constant, Yk = beta1*Xk - beta2*Xk_1.
%%
%% Input: nr = size(X,1) 
%%        nc = size(X,2)
%%        problem_type = 'NNLS' for nuclear norm regularized LS problem
%%                     = 'SDLS' for regularized semidefinite LS problem
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
             APGL(nr,nc,problem_type,Amap,ATmap,bb,mutarget,rhotarget,par)

   clear global  %% important 
   global X Xold Xh Xhold Grad 

   if ~exist('par'); par = []; end
   if ~strcmp(problem_type,'SDLS') & ~strcmp(problem_type,'NNLS')
      error([' problem_type must be ''NNLS'' or ''SDLS'' ']); 
   end
   nmin = min(nr,nc);
   randnstate = randn('state'); randstate = rand('state');
   randn('state',5); rand('state',5);

   tol     = 1e-4;
   maxiter = 100;
   maxrank = min(nmin,150); 
   verbose = 1;
   plotyes = 0;
   linesearch   = 1; 
   continuation = 1; 
   continuation_scaling = 1e-4; 
   truncation       = 1;    truncation_start = 15; 
   truncation_gap   = 5; 
 
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
   if isfield(par,'truncation_start'); 
      truncation_start = par.truncation_start; end
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
      fullsvdswitch = 300;
   else
      fullsvdswitch = max(floor(nmin/10)+100,300);
   end
%%
%% initial value of the number of singular values computed
   svp = 1; svpold = svp;  
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
   obj   = 0; objlinesearch = inf; 
   normX = 0;  
   normXold = 0; 
%%
   param.nr = nr; 
   param.nc = nc; 
   param.plotyes = plotyes; 
   param.truncation = truncation;
   param.maxrank = maxrank; 
   param.matrix_format = matrix_format;  
   param.use_fullsvd = 0; 
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
      AY_b = AY-bb; 
      Grad = ATmap(AY_b); 
      normGrad = norm(Grad,'fro'); 
      %%
      %% G = Y-Grad/tau, where Y = beta1*X-beta2*Xold. 
      %%
      if strcmp(matrix_format,'standard')
         trXXold = sum(sum(X.*Xold));     
      else
         trXXold = sum(sum((Xhold.U'*Xh.U).*(Xhold.V'*Xh.V)));
      end
      trXGrad    = AX'*(AY_b);
      trXoldGrad = AXold'*(AY_b);
      normY2  = beta1^2*normX^2 + beta2^2*normXold^2 - 2*beta1*beta2*trXXold;
      trYGrad = beta1*trXGrad - beta2*trXoldGrad;
      normG2  = normY2 - (2/tau)*trYGrad + (1/tau^2)*normGrad^2;  
      %%  
      %% preditermined number of singular values
      %%
      svold = sv; 
      param.beta1 = beta1; param.beta2 = beta2; param.beta3 = 1/tau;
      if (iter>1); param.verbose=verbose; else; param.verbose=0; end
      if (sv > fullsvdswitch); 
         param.use_fullsvd = 1; 
      end
      if (svp == sv) & (abs(mu/mutarget) < 50) & (obj < 10*objlinesearch)
         sv = min(svp+5,nmin); 
      else
         sv = min(svp+1,nmin);
         if (iter>1 & iter<15) & (svp == svpold); 
            %% delay truncation when there is no change in rank
            %% of X in the early phase of iteration.
            %% 2009-Nov-03
            truncation_start = truncation_start+1; 
            fprintf('*');
         end
      end
      svpold = svp; 
      if (abs(mu) > abs(mutarget));
         svdtol = 0.1*min([1e-3,1e-2*abs(mu)]);  
      else
         svdtol = 0.1*min([1e-5,1e-1*tol,1e-2*abs(mu)]); 
      end
      if (iter <= truncation_start) ...
         | (max(runhist.sv)<=10 & iter<=2*truncation_start) %% 2009-Nov-03
         gap = inf; 
         fprintf('+');
      else
         if (mu > 10*mutarget); %% 2012-Jan-07
            gap = 2*truncation_gap;
         else
            gap = truncation_gap;
         end
      end
      if strcmp(problem_type,'NNLS')
         [U,V,sd,svp] = proxmap(sv,1,mu/tau,svdtol,gap,param);
      else
         [U,V,sd,svp] = proxmapsym(sv,1,mu/tau,svdtol,gap,param);
      end
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
      trXnewGrad = AX'*(AY_b); 
      trXnewY = beta1*trXnewX - beta2*trXnewXold; 
      trXnewG = trXnewY - (1/tau)*trXnewGrad;      
      normXG2 = norm(sd)^2 + normG2 - 2*trXnewG; 
      normAY_b = sqrt(max(0,trYGrad - bb'*(AY_b))); 
      objlinesearch = 0.5*normAY_b^2 - normGrad^2/(2*tau) + (tau/2)*normXG2;  
      objlinesearch = objlinesearch + 0.5*rho*normY2 + mu*sum(sd); 
      %%
      if strcmp(matrix_format,'standard')
         Xold = X; X = Xnew;
      else
         Xhold = Xh; Xh = Xhnew; 
      end
      trXXold   = trXnewX; %% important
      normX     = norm(sd);
      normRes   = norm(AX-bb); 
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
      runhist.relDist(iter)       = relDist;
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
      if (mu == mutarget) & (iter > 4) & (obj < 10*objlinesearch)
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
         if (all(abs(relResratio-1) < 5*tol) | (relObjdiff < 0.1*tol)) ...
            & (relDist < 10*tol)
	    msg = sprintf('[c] lack of progress in relRes');
            fprintf('\n %s : %2.1e',msg,max(relResratio-1)); 
            break; 
         end
	 if (relResratio(end) < 0.1) & (iter > 20) %% added: 2011-Apr-01
            msg = sprintf('[d] sudden jump in relRes'); 
            fprintf('\n %s : %2.1e',msg,relResratio(end));
            if strcmp(matrix_format,'standard')
               X = Xold;   
            else
               Xh = Xhold;   
            end
            relDist = runhist.relDist(iter-1); 
            svp = svpold; 
            break;  
         end
      end
      %%
      %% perform linesearch to update tau. 
      %% 
      if (abs(mu) < 100*abs(mutarget)) & (linesearch) 
         if (verbose); fprintf(' %5.4e|',objlinesearch); end
         const_eta = 0.8; 
         if (obj < 0.9999*objlinesearch) %% old factor = 1.000; 
   	    tau = min(taumax,max(const_eta*tau,taumin)); 
	 elseif (obj > objlinesearch)
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
         if (mu == mutarget) & (iter > 4)
            idx = [iter-4:iter]; 
            if all(runhist.obj(idx) < 0.98*runhist.objlinesearch(idx))
               taumin = taumin*const_eta; 
            end
         end
      end
      %%    
      %% update mu
      %%
      mu = sign(mutarget)*max(0.7*abs(mu),abs(mutarget)); 
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

