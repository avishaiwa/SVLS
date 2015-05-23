n=[1000 1000]; r=10; k=12; noise=0.05; X_type = 'low_rank' ;measurement_type = 'columns_and_rows_normal';
long=100; SVLS_p_RMSE=zeros(long,1); SVLS_p_QR_RMSE = SVLS_p_RMSE;  SVLS_RMSE=SVLS_p_RMSE;
for i=1:long
    i
    [ m, measure ] = samp_matrix( ...
        n,r,k,noise, X_type,measurement_type  );
    
    Bc = measure.Bc; Ac = measure.Ac; Br = measure.Br; Ar = measure.Ar;
    
    N_QR= SVLS_p( Br,Bc,Ar,Ac,r, 1);
        SVLS_p_QR_RMSE(i) = norm(N_QR-m,'fro')/norm(m,'fro');
    N= SVLS_p( Br,Bc,Ar,Ac,r);
        SVLS_p_RMSE(i) = norm(N-m,'fro')/norm(m,'fro');

        
    [X] = SVLS( Br,Bc,Ar,Ac,r);
    SVLS_RMSE(i)=norm(m-X,'fro')/norm(m,'fro');
    
end

%Holk_RMSE(Holk_RMSE>1)=1;
%SVLS_RMSE(SVLS_RMSE>1)=1;

% loss = [SVLS_RMSE,Holk_RMSE];
% nbins=10;
% figure; hist(loss,nbins)
% figure; hist(Holk_RMSE-SVLS_RMSE,nbins)
 figure; plot(SVLS_p_RMSE , SVLS_RMSE,'+'); hold on; plot(0:1, 0:1, 'r'); 
 xlabel('SVLS_P RRMSE'); ylabel('SVLS RRMSE'); 
 xlim([0 1]); ylim([0 1]);
 
 
data_file_name = fullfile(['results_vec\SVLS_p' '.mat']);
save( data_file_name  , 'SVLS_p_RMSE','SVLS_RMSE');
 



