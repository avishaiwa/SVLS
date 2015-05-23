% NOTE: this should be merged with the othe simulation file (single_COLROW)

%%%symoletion for matrix completion with linear transofrmetion of all measurment. 
%%for each p and each algorithm find the RMSE ||N_hat-m||/||m|| when m is
%%nXn matrix from rank r and N_hat is the estmetion of m. 
n=40;
max_p = 400;
r = 2;
noise=0.25;
max_iter = 5000;%maximum iteretion for each algorithm
tol=10^-3;


result_normal = zeros(iter,max_p);%RMSE for gradient descent
result_APG_normal= zeros(iter,max_p);%RMSE for APG
dist_normal = zeros(iter,max_p);%number of iteretion for gradient descent
dist_APG_normal= zeros(iter,max_p);%number of iteretion for APG

for p=100:max_p
    p
    for i=1:iter
        i
        [map, tmap, am, m] = Normal_measurment(n,r,p,noise) ;
        ra = randn(n);
        [X S Y ] = OptSpace_col_func(ra,am,map,tmap,r,1000,tol);
        result_normal(i,p) = norm(X*S*Y'-m,'fro')/norm(m,'fro');
       
        [X,numiter,ttime,sd,runhist] =  APGL(40,40,'NNLS',map,tmap,am,300,10^-6,0); % minimize nuclear norm
        result_APG_normal(i,p) = norm(X-m,'fro')/norm(m,'fro');
        dist_APG_normal(iter,40) = numiter;
        
    end
    
end



