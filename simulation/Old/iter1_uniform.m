function [time,result,dist ] = iter1_uniform( m,am,alg_str,p,max_iter)
n = length(m(:,1));
switch alg_str
    case 'OptSpacenoTrim'
        tic
        [X S Y dist] = OptSpacenoTrim(am,r,max_iter,10^-2); %normal optSpace
        time =  toc;
        N_hat= X*S*Y';
        result = norm(N_hat-m,'fro')/norm(m,'fro');
        dist = length(dist);

    case 'SVT'
         tau = 5*n;
         p1 = 2*n*p -p^2;
         delta = 1.2*p1/n^2; 
          tic
          [X,S,Y,numiter,out]  = SVT(n,find(am~=0),am(am~=0),tau,delta,max_iter,10^-4);%svt algorithm
          time =  toc;
          N_hat= X*S*Y';
          result =  norm(N_hat'-m,'fro')/norm(m,'fro');
          dist = numiter;
          
    case 'OptSpace'
          tic
           [X S Y dist] = OptSpace(am,r,max_iter,10^-3); %normal optSpace
           time =  toc;
           N_hat= X*S*Y';
           result = norm(N_hat-m,'fro')/norm(m,'fro');
           dist = length(dist);
end         
                
end

