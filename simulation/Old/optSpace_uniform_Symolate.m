%function [ result] = optSpace_uniform_Symolate( iter,n,r,noise )
%input 
%iter - number of iteretion per number of measurment
%n - size of the matrix we sample
%r - rank of the matrix
%noise - m = N + z, when Z~N(0,noise)
%output
%result - matrix  with rowsre present the four algorithm for matrix
%completio and colmun p reprsent the probality tu success with p unifrom distrbuted  entries.
result_opt = zeros(iter,1600);result_optRand= zeros(iter,1600);result_grad= zeros(iter,1600);result_noS= zeros(iter,1600);

for p=874:1000
  p
    for i=1:iter
        [ am, m] = opt_measurment(n,r,p,noise);
        
        [X S Y] = OptSpace(am,r,1000,10^-4); %normal optSpace
        %if( norm(X*S*Y'-m,'fro')/norm(m,'fro')<10^-2 ) 
            result_opt(i,p) = norm(X*S*Y'-m,'fro')/norm(m,'fro');
        %end
        
        [X S Y] = OptSpaceRand(am,r,1000,10^-4);%optSpace with no S
        %if( norm(X*S*Y'-m,'fro')/norm(m,'fro')<10^-2  )
           result_optRand(i,p) = norm(X*S*Y'-m,'fro')/norm(m,'fro');
        
        %end        
        
        [X S Y] = grad_opt(randn(n),am,r,1000,10^-4);%optSpace with random start
        %if( norm(X*S*Y'-m,'fro')/norm(m,'fro')<10^-1)  
            result_grad(i,p) =norm(X*S*Y'-m,'fro')/norm(m,'fro');
        %end
        
        [X S Y] = rand_opt_noS(randn(n),am,r,1000,10^-4);%optSpace with random start and no s
        %if( norm(X*S*Y'-m,'fro')/norm(m,'fro')<10^-1)  
            result_noS(i,p) = norm(X*S*Y'-m,'fro')/norm(m,'fro');
        %end
         
    end
    %result(:,p) = result(:,p)/iter;
  
end

end

