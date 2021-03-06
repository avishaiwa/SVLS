function [ m, am,Bc,Br,Ac,Ar ] = samp_matrix(n,r,p,noise, X_type,measurement_type)
switch  measurement_type 
    case 'row_column'
        [ m,am ] =ColandRow_opt( n,r,p,noise,X_type);
        [ Br,Bc,Ar,Ac ] =optTObasis( n,p,am );
    case  'random_entries'
         [ m, am] = opt_measurment(n,r,p,noise,X_type);
    case 'ColRow_Normal_measurements'
         [m, Br, Bc, Ac, Ar] = CAR_G(n, r, p, noise, X_type);
end    
        
end

