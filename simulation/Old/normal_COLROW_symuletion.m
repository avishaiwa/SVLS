%%%symoletion for matrix completion with linear transofrmetion of col and row 
%%for each p and each algorithm find the RMSE ||N_hat-m||/||m|| when m is
%%nXn matrix from rank r and N_hat is the estmetion of m. 
n=40; % matrix size (we assume a square matrix) 
p_max = 400; % number of row and column measurements 
r = 2; % rank of unknown matrix 
noise=0.25; % st.d. of noise 
max_iter = 5000; % maximum # iterations for each algorithm
tol=10^-3;

result_COLROW_normal = zeros(iter,p_max);%RMSE for basis algorithm
dist_COLROW_normal = zeros(iter,p_max);%number of iteretion in basis algorithm

for p=100:p_max
    p
    for i=1:iter
        i
        [ m,Br,Bc,ar,ac ] = CAR_G( n,r,p,noise ); %linear transoformation of columns and rows
        [ N] = basis( Br,Bc,ar,ac,r ); % problem inside basis function !!! 
        [map,tmap,am,C] = colrow_map(ar,ac,n,p,m,Br,Bc); %change the problem to a dual problem, instert Bc=ac*m and Br= ar*m get map(m)=am.
        [X S Y dist] = OptSpace_col_func(N,am,map,tmap,r,max_iter,tol); %optSpace with the output of basis algorithm as starting point.
        result_COLROW_normal(i,p) =  norm(X*Y'-m,'fro')/norm(m,'fro');
        dist_COLROW_normal(i,p) = length(dist);       
        
    end
    
end



