% Sample a random low-rank marix, and also sample noisy measurements
% TODO: Add explanation on each variable, see if we can unite format over
% different types of measurements
%
% Input:
% n - vector sampled matrix from size n1Xn2
% r - Integer, the rank of the sampeled matrix
% k - Integer, (n1+n2)*k is the number of measurements 
% noise - Double,  add noise to the measurements
% X_type - String, low rank or symmetric low rank
% measurement_type - String:
%                       columns_and_rows- For RCMC model
%                       random_entries - standard matrix completion
%                       columns_and_rows_normal- For GRC model
%                       gaussian- For gaussian ensamble model
% Output:
% m - n1Xn2 random matrix of rank r. 
% measure.am - RCMC and standard matrix completion,  
%                          Sparse matrix with zeroes at the unrevealed indices.
% measure.Ar -  kXn2  matrix, Rows affine transformation  for RCMC and GRC model
% measure.Ac - n1XK  matrix, Columns affine transformation  for RCMC and GRC model
% measure.Br - Ar*m+N(0,noise^2), Rows measurements for RCMC and GRC model
% measure.Bc - m*Ac+N(0,noise^2), Columns measurements for RCMC and GRC model
% measure.map - function handle, linear transformation for gaussian ensable model
% measure.tmap - function handle, transpose of map (why needed?) 
% measure.bb -  vector of length (n1+n2)k, bb=map(m) + N(0,noise^2), (noisy version of map)
% measure.n1, measure.n2, measure.k, measure.k1 - parameers for sampling 
%
function [ m, measure ] = samp_matrix( ...
    n, r, k, noise, X_type, measurement_type )
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
        m = u*v' + 0*randn(n1,n2);
    case 'symmetric_low_rank'
        lambda = diag(randn(r,1));
        lambda = lambda/norm(lambda,'fro');
        m = u*lambda*u';
    case 'positive_definite_low_rank'
        m=u*u'; 
    
end

% Next, sample measurements
switch lower(measurement_type)
    case 'columns_and_rows'
        k1 = (n1+n2)*k-k^2; % compute total number of scalar measurements
        measure.am =  ColandRow_sample(m, k, noise, X_type); % X_type isn't used here!
        [measure.Br, measure.Bc, measure.Ar, measure.Ac ] = create_measurements(n, k, measure.am); % What's here? why sample again???
    case  'random_entries'
        k1 = (n1+n2)*k-k^2; % compute total number of scalar measurements
        measure.am = random_measurment(m, k1, noise, X_type);
    case {'columns_and_rows_normal', 'gaussian_columns_and_rows'}
        k1 = (n1+n2)*k-k^2; % compute total number of scalar measurements
        [measure.Br, measure.Bc, measure.Ac, measure.Ar] = CAR_GRC(m, k, noise, X_type);
    case 'gaussian'
        k1 = (n1+n2)*k; % compute total number of scalar measurements (why different from 'random_entries'?)
        [measure.map, measure.tmap, measure.bb] = Normal_measurment(m, n,k1, noise); % why different variables?
        measure.am = measure.bb; % rename field?? 
    case 'rank_one_projection' % method of Cai's paper (save map? save only vector gamma and beta)
        % here Ar, Ac are collections of vectors for ROP measurements.
        % Different meaning from column-and-row
        [measure.am, measure.Ar, measure.Ac] = ROP_measurement(m, k, noise, X_type);
end

measure.n1 = n1; measure.n2 = n2; measure.k = k; measure.k1 = k1; % record also parameters

