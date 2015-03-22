%%************************************************************************
%% run random matrix completion problems. 
%% 
%% 
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
%%************************************************************************

  clear; 

  addpath('PROPACKmod/');
  addpath('solver');
%%
%% generate random a test problem 
%% n = size of the unknown square matrix M, r = rank of M
%%
  ndim  = [1e3, 1e3, 1e3, 5e3, 5e3, 5e3, 1e4, 1e4, 1e4, 2e4, 3e4, 5e4, 1e5];
  rrank = [10,  50,  100, 10,  50,  100, 10,  50,  100, 10,  10,  10,  10]; 
  pfac  = [6,   4,   3,   6,   5,   4,   6,   5,   4,   6,   6,   6,   6];
%%
  symmetric = 0; 
  %%scenario = 'noiselss';  
  scenario = 'noisy'; 
  randstate_idx = [1:2]; 
  par.tol     = 1e-4;
  par.verbose = 1;
  par.plotyes = 0; 
%%
  for kk = [1]; 
      r  = rrank(kk); 
      cp = pfac(kk); 
      n  = ndim(kk);
      nr = n; nc = n;  
      if strcmp(scenario,'noiseless')
         noiseratio = 0;   %% noiseless case     
      else
         noiseratio = 0.1; %% noisy case
      end
      for randstate = [randstate_idx];   
         fprintf('\n Create matrix %2.0d with rank = %2.0d,',randstate,r); 
         randn('state',double(randstate));
         rand('state',double(randstate));  
         if (symmetric)
            dr = n*r - r*(r-1)/2; 
             p  = 2*cp*dr; %% number of sampled entries    
         else
            dr = r*(nr+nc-r);
            p  = cp*dr; %% number of sampled entries    
         end
         %% construct M (nxn), Omega (nxn)
         %% M = ML*MR';
         %% Omega = spconvert([II,JJ,ones(p,1); n,n,0]); 
         %% B = M.*Omega; 
         ML = randn(nr,r); 
         if (symmetric)
            MR = ML;
         else
            MR = randn(nc,r); 
         end
         normM = sqrt(sum(sum((ML'*ML).*(MR'*MR))));
         prob = p/(nr*nc); 
         II = zeros(p,1); JJ = zeros(p,1); cnt = 0; 
         if (symmetric)
            for j=1:nc
               tmp = rand(j-1,1); 
               idx = find(tmp < prob); 
               idxkeep = [idx; j];   
               len = length(idxkeep); 
               II(cnt+[1:len]) = idxkeep;  
               JJ(cnt+[1:len]) = j*ones(len,1);  
               cnt = cnt + len; 
            end
         else
            for j=1:nc
               tmp = rand(nr,1); 
               idx = find(tmp < prob);  
               len = length(idx); 
               II(cnt+[1:len]) = idx;  
               JJ(cnt+[1:len]) = j*ones(len,1);  
               cnt = cnt + len; 
            end
         end
         II = II(1:cnt); JJ = JJ(1:cnt); p = cnt; 
         %%
         Jcol = compJcol(JJ); 
         bb = ProjOmega(ML,MR,II,Jcol); 
         B  = spconvert([II,JJ,bb; n,n,0]); 
         if strcmp(scenario,'noiseless')
            xi = sparse(p,1); 
            sigma = 0; 
         else
            randnvec = randn(p,1);
            sigma = noiseratio*norm(bb)/norm(randnvec); 
            xi = sigma*randnvec; 
            B  = B + spconvert([II,JJ,xi; nr,nc,0]); 
         end
         fprintf(' noise ratio: sigma = %3.2e\n',sigma); 
         %%------------------------------------------------
         %% evaluate the regularization parameter mu
         %%
         options.tol = 1e-8; 
         mumax = svds(sparse(B),1,'L',options);
         const = 1e-4;  
         mutarget = const*mumax;
         noiseratio = norm(xi)/norm(bb); 
         %%
         fprintf('\n Problem: n = %d, p = %d, r = %d,',n,p,r); 
         fprintf(' p/dr = %3.2e, p/n^2 = %3.2e%%',p/dr, p/(n^2)*100)
         fprintf('\n mu = %3.2e, noise ratio: sigma = %3.2e,',mutarget,sigma); 
         tstart = clock;
         if (symmetric)
            B = B+triu(B,1)'; 
            [X,iter,time,sd,hist] = APGLsym(B,mutarget,0,par);
         else
            [X,iter,time,sd,hist]= APGL(B,mutarget,0,par);
         end
         if isstruct(X)
            normX = sqrt(sum(sum((X.U'*X.U).*(X.V'*X.V))));
            trXM = sum(sum((ML'*X.U).*(MR'*X.V)));
         else
            normX = norm(X,'fro'); trXM = sum(sum(ML.*(X*MR))); 
         end
         mse = sqrt(normX^2+normM^2-2*trXM)/normM;
         runhist.mu(randstate)     = mutarget; 
         runhist.mumax(randstate)  = mumax; 
         runhist.time(randstate)   = etime(clock,tstart);
         runhist.iter(randstate)   = iter;
         runhist.obj(randstate)    = hist.obj(end);
         runhist.mse(randstate)    = mse; 
         runhist.svp(randstate)    = hist.svp(end);
         runhist.maxsvp(randstate) = max(hist.svp); 
         runhist.maxsig(randstate) = max(sd); 
         runhist.minsig(randstate) = min(sd(find(sd>0)));
         %%
         %% report results in a table
         %%
         fprintf(' Problem: nr = %d, nc = %d, p = %d, r = %d,',nr,nc,p,r); 
         fprintf(' p/dr = %3.2e, p/n^2 = %3.2e%%',p/dr, p/(n^2)*100)
         fprintf('\n mu = %3.2e, noise ratio: sigma = %3.2e',mutarget,sigma)
         fprintf('\n-----------------------------------------------');
         fprintf('------------------------------')
         fprintf('\n randstate          :  %6.0d',randstate);
         fprintf('\n iterations         :  %6.0d',runhist.iter(randstate));
         fprintf('\n # singular         :  %6.0d',runhist.svp(randstate));
         fprintf('\n obj  value         :  %6.5e',runhist.obj(randstate));
         fprintf('\n cpu   time         :  %6.2e',runhist.time(randstate));
         fprintf('\n mean square error  :  %6.2e',runhist.mse(randstate));
         fprintf('\n------------------------------------------------'); 
         fprintf('------------------------------\n')
      end
  end
%%*************************************************************************
