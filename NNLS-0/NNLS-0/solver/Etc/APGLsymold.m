%%**********************************************************************
%% Feb, 2009
%%  This is a Matlab implementation of an accelerated proximal gradient 
%%  singular value thresholding algorithm for solving the 
%%% matrix completion problem
%%       min_X  f(X) + mu||X||_*,
%% where f(X) = 0.5||A(X)-b||^2_2 + 0.5*rho||X||^2_F.
%%
%% At each iteration k, we solve the following subproblem:
%%       min_X  0.5*L||X-(Y^k-\nabla f(Y^k)/L)||^2_F + mu||X||_*,
%% where L is a Lipschitz constant, Y^k = X^k + t^k(1/t^{k-1}-1)(X^k-X^{k-1}),
%% and 0<t^k<=1, to get X^k.
%%**********************************************************************

   function [X,iter,ttime,sd,runhist] = ...
             APGLsym(const,B,mutarget,rhotarget,par)

   clear global  %% important 
   global Xh Xhold Grad beta1 beta2 tau 

   if ~exist('par'); par = []; end

   [nr,nc] = size(B);
   if (norm(B-B','fro'))
      error(' B must be a symmetric matrix'); 
   end
   nmin = min(nr,nc);

   tol     = 1e-3;
   maxiter = 100;
   maxrank = min(nmin,200); 
   verbose = 1;
   plotyes = 0;
   linesearch   = 1; 
   continuation = 1; 
   truncation   = 1; 
 
   tstart = clock;
   if isfield(par,'tol'); tol = par.tol; end
   if isfield(par,'maxiter'); maxiter = par.maxiter; end
   if isfield(par,'maxrank'); maxrank = par.maxrank; end
   if isfield(par,'verbose'); verbose = par.verbose; end
   if isfield(par,'plotyes'); plotyes = par.plotyes; end
   if isfield(par,'linesearch');   linesearch = par.linesearch; end   
   if isfield(par,'continuation'); continuation = par.continuation; end
   if isfield(par,'truncation');   truncation = par.truncation; end

%% initial mu and target mu
   if (~continuation)
      const = 1; 
   end
   mu  = mutarget/const;
   rho = rhotarget/const; 
%%
   if (nmin <= 1000)
      fullsvdtol = 200;
   else
      fullsvdtol = max(floor(nmin/10)+100,200);
   end
%%
%% initial value of the number of singular values computed
   svp = 1; 
   sv  = svp; 
%%
%% Initialization
%%
   if (max(nr,nc) < 1000);  
      SVDoptions = 1; 
   else
      SVDoptions = 2; 
   end   
   if (verbose); 
      fprintf('\n\n  nr = %2.0d, nc = %2.0d, nnz = %2.0d',nr,nc,nnz(B)); 
      fprintf('\n  mu_target = %3.2e',mutarget); 
      fprintf('\n  SVD options = %2.0d',SVDoptions); 
      fprintf('\n-----------------------------------------------');
      fprintf('-------------------------------------------'); 
      fprintf('\n it     mu       tau      svp    relXdiff  relRes '); 
      fprintf(' relObjdiff   obj       time')
      fprintf('\n-----------------------------------------------');
      fprintf('-------------------------------------------');  
   end   
   [II,JJ,bb] = find(B);
   Jcol = compJcol(JJ); 
   clear B;  
   if (SVDoptions == 1)
      X = sparse(nr,nc);
      Xold = X; 
   else
      Xh.U = sparse(nr,1); Xh.V = sparse(nc,1);
      Xhold = Xh;
   end
   AX = sparse(length(bb),1); 
   AXold = AX;
   normb = max(norm(bb,'fro'),1); 
   obj   = 0;
   normX = 0;  
   normXold = 0; 
%% 
   taumax = 1; tau = taumax; taumin = 1e-3; 
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
         AY    = (1+beta)*AX - beta*AXold; 
      end
      Gradvec = AY-bb; 
      Grad = spconvert([II,JJ,Gradvec; nr,nc,0]); 
      normGrad = norm(Grad,'fro'); 
      %%
      %% G = Y-Grad/tau, where Y = beta1*X-beta2*Xold. 
      %%
      if (SVDoptions == 1)
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
      if (svp == sv) & (mu/mutarget < 50)
         sv = min(svp+5,nmin); 
      else
         sv = min(svp+1,nmin);
      end
      if (SVDoptions == 1); HH = beta1*X - beta2*Xold - Grad/tau; end
      if (sv > fullsvdtol) & (SVDoptions == 1)
         countmax = 1; 
      end
      countmax = 1;
      count = 1; 
      while (count<=countmax) 
         if (mu > mutarget);
            svdtol = min(1e-2,0.01*mu);  
         else
            svdtol = min(0.1*tol,0.01*mu); 
         end
         options.tol = svdtol; 
         if (SVDoptions == 1)
            if (sv > fullsvdtol) 
               [U,S] = mexeig(full(HH));
               countmax = 1; 
            else
               [U,S] = laneig(HH,sv,'AL',options);
            end
         else
            [U,S] = laneig('Amatvec',nr,sv,'AL',options);
         end
         d = diag(S);
         sd = max(0,d-mu/tau);
         %%
         %% estimate the rank 
         %%
         if (iter <= 10); 
            gap = inf; 
         else
            gap = 5; 
         end
         ratio  = d(1:end-1)./d(2:end);
         idxbig = find(ratio > gap);   
         if ~isempty(idxbig)
            idxtruncate = idxbig(1)+1;          
         else
            idxtruncate = length(find(sd>0)); 
         end
         normratio = zeros(1,idxtruncate); 
         window = 10; 
         for kk = 1:idxtruncate-1
            sd1 = sd(max(1,kk-window):kk); 
            sd2 = sd(kk+1:min(idxtruncate,kk+window+1)); 
            normratio(kk) = mean(sd1)/(mean(sd2)+eps); 
         end 
         [maxnormratio,maxidx] = max(normratio); 
         if (verbose & iter > 1); fprintf(' %2.1e',maxnormratio); end
         if (truncation)
            if (maxnormratio > gap)
               sd = sd(1:maxidx); 
            end
         end
         sd = sd(1:min(length(sd),maxrank)); 
         %% 
         %% plot graph
         %%
         if (plotyes)
            subplot(121)
            semilogy([1,length(d)],mu/tau*[1,1]); hold on;               
            semilogy([1,length(d)],1e1*mu/tau*[1,1]);
            semilogy([1,length(d)],1e2*mu/tau*[1,1]);
            semilogy(d,'*'); hold off;  
            subplot(122)
            plot(ratio,'*'); pause(0.01);
         end
         %%
         %% there may be more singular values that are >= mu/tau.
         %% So increase the predetermined number of singular values.
         %%
         svp = length(find(sd>0));
         if (svp < sv) | (svp == nmin); 
            break; 
         else
            count = count + 1; 
         end
      end 
      Shlf = spdiags(sqrt(sd),0,svp,svp); 
      U = full(U(:,1:svp)*Shlf); 
      V = U; 
      %%
      if (SVDoptions == 1)
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
      if (exist('mexProjOmega') == 3)
        AX = mexProjOmega(U,V,II,Jcol); %% Note: must be after AXold = AX; 
      else
        AX = ProjOmega(U,V,II,JJ); %% Note: must be after AXold = AX; 
      end
      trXnewGrad = AX'*Gradvec; 
      trXnewG = beta1*trXnewX - beta2*trXnewXold - (1/tau)*trXnewGrad;      
      normXG2 = norm(sd)^2 + normG2 - 2*trXnewG; 
      objlinesearch = 0.5*normGrad^2 - normGrad^2/(2*tau) + (tau/2)*normXG2;  
      objlinesearch = objlinesearch + 0.5*rho*normY2 + mu*sum(sd); 
      %%
      if (SVDoptions == 1)
         Xold  = X; 
         X     = U*V';  
         trXXold = sum(sum(X.*Xold));
      else
         Xhold = Xh; 
	 Xh.U  = U; Xh.V = V; 
         trXXold = sum(sum((Xhold.U'*U).*(Xhold.V'*V)));
      end      
      normX       = norm(sd,'fro');
      normRes     = norm(AX-bb,'fro'); 
      obj         = 0.5*normRes^2 + mu*sum(sd) + 0.5*rho*normX^2;
      normXdiff   = sqrt(normX^2 + normXold^2 - 2*trXXold); 
      relXdiff   = normXdiff/max(normX,1);
      relObjdiff = abs(obj-objold)/max(1,obj);  
      relRes     = normRes/normb; 
      runhist.obj(iter)           = obj; 
      runhist.objlinesearch(iter) = objlinesearch; 
      runhist.relXdiff(iter)      = relXdiff;
      runhist.relObjdiff(iter)    = relObjdiff;  
      runhist.relRes(iter)        = relRes;  
      runhist.tau(iter)           = tau;
      runhist.mu(iter)            = mu; 
      runhist.svp(iter)           = svp; 
      runhist.sv(iter)            = sv; 
      %%
      if (verbose)
         fprintf('\n %2.0d  %3.2e %3.2e |%3.0d %3.0d| ',iter,mu,tau,svp,sv);
	 fprintf('%3.2e %3.2e %3.2e| %5.4e|',relXdiff,relRes,relObjdiff,obj); 
      end    	  
      %% 
      %% check stopping criterion
      %%
      if (mu == mutarget) & (iter > 4)
         if (max([relXdiff, 0.1*relObjdiff]) < tol) 
	    msg = sprintf('relXdiff < %3.2e',tol);
 	    fprintf('\n %s',msg); 
            break; 
         end
	 if (max([0.5*relXdiff, 100*relObjdiff]) < tol)  
	    msg = sprintf('relObjdiff < %3.2e',0.01*tol);
 	    fprintf('\n %s',msg); 
            break; 
         end
         idx = [iter-2:iter];
         relResratio = runhist.relRes(idx-1)./runhist.relRes(idx); 
         if all(abs(relResratio-1) < 5*tol) ...
            & (relObjdiff < 0.1*tol) & (relXdiff < 10*tol)
	    msg = sprintf('lack of progress in relRes');
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
            if (SVDoptions == 1)
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
      mu = max(0.7*mu,mutarget); 
      if (verbose); fprintf(' %3.2e|',etime(clock,tstart)); end
   end 
   sinmax = max(sd);
   sinmin = min(sd(find(sd>0)));
   ttime  = etime(clock,tstart);
   if (SVDoptions==2)
      X = Xh; 
   end
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
      fprintf(1,' norm(X-Xold,''fro'')/norm(X,''fro'') = %3.2e\n',relXdiff);
      fprintf(1,'\n');
   end
%%**********************************************************************

