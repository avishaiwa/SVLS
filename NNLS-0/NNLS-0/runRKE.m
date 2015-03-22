%%************************************************************************
%% Kernal regularization problem: 
%% 
%% min_{X psd} sum_{ij\in nnz(W)} (<Aij,X>-dij)^2 + mu*Tr(X)
%% where Aij = eij*eij', eij = ei-ej.
%%
%% For details, see: 
%% [1] Lu, F., Keles, S., Wright, S., and Wahba, G. 
%%     Framework for Kernel Regularization With Application to Protein 
%%     Clustering, 
%%     Proceedings of the National Academy of Sciences, 102 (2005).
%%
%% NNLS, version 0:
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun
%%************************************************************************

   clear all

   addpath('RKEData'); 
%%
   load proteinDisMatrix
   n = 630; 
   dtmp = d(1:n,1:n); 
   numneighbor = n; %% old: min(50,n); 
   if (numneighbor == n)
      H = ones(n,n); 
   else
      H = zeros(n,n); 
      for j=1:n
         perm = randperm(n); 
         idx = perm(1:numneighbor); 
         H(idx,j) = ones(length(idx),1); 
      end
      H = min(1,H + H'); 
   end   
%%
   tic; 
   [blk,At,b,numeq] = RKESDPdata(dtmp,H); 
   toc; 
   randstate = randn('state'); 
   randn('state',0);
   eigsopts.disp  = 0; 
   eigsopts.tol   = 1e-6; 
   eigsopts.maxit = 200;
   eigsopts.v0    = randn(n,1); 
   randn('state',randstate);
   Atb = mexsmat(blk,At{1}*b,1);
   [V,D,flag] = eigs(Atb,3,'LM',eigsopts); 
   mumax    = max(abs(diag(D))); 
   mufactor = 1e-3; 
   mutarget = mufactor*mumax; 
%%
   Amap  = @(X) (mexsvec(blk,X)'*At{1})'; 
   ATmap = @(y) mexsmat(blk,At{1}*y,1); 
   par.tol     = 1e-4;    
   par.plotyes = 1; 
   par.truncation       = 1; 
   par.truncation_gap   = 5.0; 
   par.continuation_scaling = mufactor; 
   problem_type = 'SDLS'; 
   figure(1)
   [X,iter,ttime,sd,runhist] = ...
      APGL(n,n,problem_type,Amap,ATmap,b,mutarget,0,par);
   [V,sig] = mexeig(X); sig = full(diag(sig)); 
   sig = abs(sig); sig = sig(end:-1:1); V = V(:,end:-1:1); 
   fprintf('\n RKE: n = %3.0f, mutarget = %3.2e',n,mutarget); 
   fprintf('\n trace(X) = %3.2e\n',sum(sig));
%%
   set(0,'defaultaxesfontsize',18);
   set(0,'defaultlinemarkersize',10);
   set(0,'defaulttextfontsize',18);
   X3 = V(:,1:3)*diag(sqrt(sig(1:3)));
   figure(2) 
   plot3(X3(:,1),X3(:,2),X3(:,3),'.');
   axis('square'); grid;
   xlabel('x'); ylabel('y'); zlabel('z')
   axis(0.5*[-1, 1, -1, 1, -1, 1])
   %%print -depsc RKE-630.eps 
%%************************************************************************
