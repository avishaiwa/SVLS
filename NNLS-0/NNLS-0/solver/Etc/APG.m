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

   function [X,obj,iter,svn,ttime] = ...
             APG(const,nr,nc,B,mutarget,rhotarget,par)

   clear global  %% important 
   global Xh Xhold Grad beta1 beta2 tau 

   if ~exist('par'); par = []; end

   tol = 1e-3;
   maxiter = 200;
   verbose = 1;
   plotyes = 1;
 
   t0 = clock;
   if isfield(par,'Tol'); tol = par.Tol; end
   if isfield(par,'Maxiter'); maxiter = par.Maxiter; end
   if isfield(par,'Verbose'); verbose = par.Verbose; end

%% initial mu and target mu
   mu  = mutarget/const;
   rho = rhotarget/const; 
   fprintf(' mu_target = %3.2e\n',mutarget')
%%
   nmin = min(nr,nc);
   if (nmin <= 1000)
      ntol = 300;
   else
      ntol = max(floor(nmin/10)+200,300);
   end
%%
%% initial value of the number of singular values computed
   if (nmin >= 100)
      svn = 10; %% old: 100
   else
      svn = nmin;
   end
%%
%% Initialization
%%
   SVDoptions = 2; 
   [II,JJ,bb] = find(B); 
   if (SVDoptions == 1)
      X = sparse(nr,nc);
      Xold = X; 
      Xmask    = sparse(nr,nc); 
      Xoldmask = sparse(nr,nc);
      mask = spconvert([II,JJ,ones(length(II),1);nr,nc,0]);  
   else
      B = bb; 
      Xh.U = sparse(nr,1); Xh.V = sparse(nc,1);
      Xhold = Xh;
      Xmask = sparse(length(bb),1); 
      Xmaskold = Xmask;
   end
   normB = norm(B,'fro'); 
   obj   = 0;
   normX = 0;  
%% 
   for iter=1:maxiter

      tau = 1; 
      if (iter == 1)
         t     = 1;
         beta1 = 1; 
         beta2 = 0; 
         Ymask = Xmask; 
      else
         t     = (1+sqrt(1+4*told^2))/2; 
         beta  = (told-1)/t; 
         beta1 = (1+beta)*(1-rho/tau); 
         beta2 = beta*(1-rho/tau); 
         Ymask = (1+beta)*Xmask - beta*Xoldmask; 
      end
      if (SVDoptions == 1)
         Grad = Ymask - B;  
         HH = beta1*X - beta2*Xold - Grad/tau;
      else
         Grad = spconvert([II,JJ,Ymask-B;nr,nc,0]); 
      end
      %%  
      %% preditermined number of singular values
      sv = min(svn+1,nmin);
      if (sv > ntol) & (SVDoptions == 1)
         [U,S,V] = svd(full(HH),0);
         d   = diag(S);
         sd  = max(0,d - mu/tau);
         svn = length(find(sd>0));
         Shlf = spdiags(sqrt(sd),0,nmin,nmin); 
         U = U*Shlf; V = V*Shlf; 
      else
         count = 1; 
         if (mu/mutarget > 100); countmax = 1; else; countmax = 2; end
         while (count<=countmax) & (sv<=nmin)          
            options.tol = min(1e-3,0.1*mu);  
            if (SVDoptions == 1)
               [U,S,V] = lansvd(HH,sv,'L',options);
            else
               [U,S,V] = lansvd('Amatvec','ATmatvec',nr,nc,sv,'L',options);
            end
            d = diag(S);
            sd = max(0,d-mu/tau);
            ratio = d(1:end-1)./d(2:end);
            [maxratio,maxidx] = max(ratio);
            fprintf(' %2.1e',maxratio) 
            if (mu/mutarget > 100); 
               gap = 3; 
            elseif (mu/mutarget > 10) 
               gap = 5; 
            else 
               gap = 10;
            end
            if (maxratio > gap) 
               sd  = sd(1:maxidx); 
            end
            svn = length(find(sd>0));
            if (plotyes)
               subplot(121)
               semilogy([1,length(d)],mu/tau*[1,1]); hold on;               
               if isempty(maxratio); const2 = 10; else; const2 = max(2*maxratio,10); end
               semilogy([1,length(d)],const2*mu/tau*[1,1]);
               semilogy(d,'*'); hold off;  
               subplot(122)
               plot(ratio,'*'); pause(0.01);
            end
            %%
            %% there may be more singular values that are >= mu/tau.
            %% So increase the predetermined number of singular values.
            %%
            if (svn < sv) | (svn == nmin); 
               break; 
            else
               count = count + 1; 
               sv = min(max(sv+5,ceil(1.1*sv)),nmin);
            end
         end 
         Shlf = spdiags(sqrt(sd),0,svn,svn); 
         U = U(:,1:svn)*Shlf; V = V(:,1:svn)*Shlf; 
      end              
      %%
      told  = t;      
      Xoldmask = Xmask; 
      normXold = normX; 
      objold   = obj;       
      if (SVDoptions==1)
         Xold  = X; 
         X     = U*V';  
         Xmask = X.*mask; 
         trXXold = sum(sum(X.*Xold));
      else
         Xhold = Xh; 
	 Xh.U  = U; Xh.V = V; 
         Xmask = maskfun(U,V,II,JJ);
         trXXold = sum(sum((Xhold.U'*U).*(Xhold.V'*V)));
      end
      normRes = norm(Xmask-B,'fro'); 
      obj   = 0.5*normRes^2 + mu*sum(sd) + 0.5*rho*normX^2;
      normX = norm(sd,'fro');
      normXdiff   = sqrt(normX^2 + normXold^2 - 2*trXXold); 
      rel_Xdiff   = normXdiff/max(normX,1);
      rel_objdiff = abs(obj-objold)/max(1,obj);  
      rel_res     = normRes/normB;      
      %%
      if (verbose)
         fprintf('\n %2.0d  %3.2e %3.2e |%3.0d %3.0d| ',iter,mu,tau,svn,sv);
	 fprintf('%3.2e %3.2e %3.2e| %5.4e|',rel_Xdiff,rel_res,rel_objdiff,obj); 
      end    	  
      %% 
      %% check stopping criterion
      %%
      if (mu <= mutarget)
         if (max([rel_Xdiff, 0.2*rel_res, 0.1*rel_objdiff]) < tol) 
	    msg = sprintf('rel_Xdiff < %3.2e',tol);
 	    fprintf('\n %s',msg); 
            break; 
         end
	 if (max([0.5*rel_Xdiff, 0.1*rel_res, 100*rel_objdiff]) < tol)  
	    msg = sprintf('rel_objdiff < %3.2e',0.01*tol);
 	    fprintf('\n %s',msg); 
            break; 
         end
      end
      %%    
      %% update mu
      %%
      if (iter<5);      rp = 0.7; 
      elseif (iter<10); rp = 0.7; 
      elseif (iter<20); rp = 0.6; 
      else;             rp = 0.6; 
      end 
      if ((mod(iter,3)==0) | (rel_Xdiff<1e-5)) & (mu>mutarget)
         mu = max(rp*mu,mutarget);
	 rho = rp*rho; 
      end       
   end 
   sinmax = max(sd);
   sinmin = min(sd(find(sd>0)));
   ttime  = etime(clock,t0);
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
      fprintf(1,' Number of singular values = %2.0d\n',svn);
      fprintf(1,' max singular value        = %3.2e\n',sinmax);
      fprintf(1,' min singular value        = %3.2e\n',sinmin);
      fprintf(1,' CPU time                  = %3.2e\n',ttime);
      fprintf(1,' norm(X-Xold,''fro'')/norm(X,''fro'') = %3.2e\n',rel_Xdiff);
      fprintf(1,'\n');
   end
%%**********************************************************************

