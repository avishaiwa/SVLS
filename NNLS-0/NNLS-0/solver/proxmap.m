%%******************************************************************
%% proxmap: solve
%%          min { 0.5*norm(X-H,'fro')^2 + rho*||X||_* }
%%
%% H = beta1*X - beta2*Xold - beta3*Grad.
%%
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun
%%******************************************************************

   function [U,V,sd,svp] = ...
             proxmap(sv,countmax,rho,svdtol,gap,param);

   global X Xold Grad

   nr = param.nr; 
   nc = param.nc; 
   nmin = min(nr,nc); 
   verbose       = param.verbose; 
   plotyes       = param.plotyes; 
   truncation    = param.truncation;
   maxrank       = param.maxrank; 
   matrix_format = param.matrix_format;  
   use_fullsvd   = param.use_fullsvd; 
%%
   if strcmp(matrix_format,'standard'); 
      HH = param.beta1*X - param.beta2*Xold - param.beta3*Grad; 
   end
%%
   count = 1; 
   while (count<=countmax) 
       options.tol = svdtol; 
       if strcmp(matrix_format,'standard')
          if (use_fullsvd) 
             [U,S,V] = mexsvd(full(HH),0);
             %% guard against erroreous decomposition produced by mexsvd
             randvec = randn(size(HH,2),1); 
             tmp = (randvec'*V)';
             err = norm(HH*randvec-U*(S*tmp))/(norm(randvec)*max(1,S(1,1)));
	     if (err > 1e-12) 
		fprintf('\n mexsvd''s result is unreliable, switch to svd'); 
                [U,S,V] = svd(full(HH),0);
             end
             countmax = 1; 
          else
             [U,S,V] = lansvd(HH,sv,'L',options);
          end
       else
          [U,S,V] = lansvd('matvec_PP','matvectransp_PP',nr,nc,sv,'L',options,param);
       end
       d = diag(S);
       sd = max(0,d-rho);
       %%
       %% estimate the rank 
       %%
       d2 = d-rho; 
       ratio  = d2(1:end-1)./d2(2:end);
       idxstart = 5; %% 2009 Oct 30
       idxbig = find(abs(ratio(idxstart:end)) > gap);   
       if ~isempty(idxbig)
          idxtruncate = (idxstart-1)+idxbig(1)+1;
       else
          idxtruncate = length(find(sd>0)); 
       end
       normratio = zeros(1,idxtruncate); 
       window = 10;
       for kk = idxstart:idxtruncate-1
          sd1 = d2(max(1,kk-window):kk); 
          sd2 = d2(kk+1:min(idxtruncate,kk+window+1)); 
          normratio(kk) = mean(sd1)/(mean(sd2)+eps); 
       end 
       [maxnormratio,maxidx] = max(abs(normratio)); 
       if (verbose); fprintf(' [%2.1e]',maxnormratio); end
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
          semilogy([1,length(d)],rho*[1,1]); hold on; 
          semilogy([1,length(d)],1e1*rho*[1,1]);
          semilogy([1,length(d)],1e2*rho*[1,1]);
          semilogy(d,'*'); hold off;  
          subplot(122)
          idxpos = find(ratio > 0); idxneg = find(ratio < 0); 
          semilogy(idxpos,ratio(idxpos),'*'); hold on; 
          semilogy(idxneg,abs(ratio(idxneg)),'or'); hold off; 
          pause(0.01);
       end
       %%
       %% there may be more singular values that are >= rho.
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
    U = full(U(:,1:svp)*Shlf); V = full(V(:,1:svp)*Shlf); 
%%******************************************************************
