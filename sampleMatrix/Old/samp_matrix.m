% Sample a random low-rank marix, and also sample noisy measurements
% TODO: Add explanation on each variable, see if we can unite format over
% different types of measurements
%
% Input: 
% n -Integer samle matrix from size xXx
% r -Integer, the rank of the sampeled matrix
% k -Integer, number of measurements
% noise -Double,  add noise to the measurements
% X_type - String, low rank or symmetric low rank
% measurement_type - String: 
%                       columns_and_rows- For RCMC model
%                       random_entries- standard matrix completion
%                       columns_and_rows_normal- For GRC model
%                       gaussian- For gaussian ensamble model
% Output: 
% m - random matrix of rank r. Each entry is .. 
% am - affine measurement matrix .. 
% Br - measurements of the rows 
% Bc 
% Ac 
% Ar 
% map
% tmap
%
function [ m, am, Br, Bc, Ac, Ar, map, tmap ] = samp_matrix( ...
n,r,k,noise, X_type,measurement_type  )
Br=0; Bc=0;Ac=0;Ar=0;m=0;am=0;map=0;tmap=0;
if length(n) == 1,
    n1 = n(1); n2 = n1;
elseif length(n) == 2,
    n1 = n(1); n2 = n(2);
end


% First sample 'true' matrix m: what is the variance of each entry? should say!
u = randn(n1,r)/r^0.25;  % Why normalization by r^0.25? 
switch X_type
    case 'low_rank'
        v = randn(n2,r)/r^0.25;
        m = u*v';
    case 'symmetric_low_rank'
        m = u*u';
end



% Next, sample measurements 
switch lower(measurement_type)
    case 'columns_and_rows'
        [ am ] =  ColandRow_sample(m, k, noise, X_type); % X_type isn't used here!
        [ Br,Bc,Ar,Ac ] = create_measurements(n, k, am); % What's here? why sample again??? 
    case  'random_entries'
        k1 = (n1+n2)*k-k^2;
        [am] = random_measurment(m, k1, noise, X_type);
    case 'columns_and_rows_normal'
        am = zeros(n);
        [Br, Bc, Ac, Ar] = CAR_GRC(m, k, noise, X_type);        
    case 'gaussian'
        k1 = (n1+n2)*k;
        [map, tmap, am] = Normal_measurment(m, n,k1, noise); % why different variables? 
        
    case 'rank_one_projection' % method of Cai's paper (save map? save only vector gamma and beta) 
        % here Ar, Ac are collections of vectors for ROP measurements. 
        % Different meaning from column-and-row 
        [am, Ar, Ac] = ROP_measurement(m, k, noise, X_type); 
        
end


