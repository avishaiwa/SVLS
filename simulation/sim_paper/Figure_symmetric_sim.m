n=500; r=5; tau=[0.1 0.25];alg_str={'SVLS_sym'};num_sampled_matrices=100; measurement_type='columns_and_rows'; k_vec = [r:25];
maximum_iter=100;tol=10^-4;
losses = ones(4,max(k_vec));
for noise=tau
    X_type='symmetric_low_rank';
    j=r;
    for k=k_vec
        k
        temp = zeros(num_sampled_matrices,1);
        for i=1:num_sampled_matrices
            [ m, measure ] = samp_matrix( ...
                n,r,k,noise, X_type,measurement_type  );
            
            [X_hat] = SVLS_sym( measure.Br, measure.Bc, measure.Ar, measure.Ac,[]);
            temp(i) = RRMSE(X_hat,m);
            
        end
        temp(temp>1)=1;
        
        mean(temp)
        if noise==0.1
            losses(1,k) = mean(temp);
        else
            losses(3,k) = mean(temp);
        end
    end
    
    X_type='low_rank';
    j=r;
    for k=k_vec
        k
        noise
        temp = zeros(num_sampled_matrices,1);
        for i=1:num_sampled_matrices
            [ m, measure ] = samp_matrix( ...
                n,r,k,noise, X_type,measurement_type  );
            
            [X_hat] = SVLS_sym( measure.Br, measure.Bc, measure.Ar, measure.Ac,[]);
            temp(i) = RRMSE(X_hat,m);
        end
        temp(temp>1)=1;
        mean(temp)
        
        if noise==0.1
            losses(2,k) = mean(temp);
        else
            losses(4,k) = mean(temp);
        end
    end
    
end
figure;  hold on;
plot(losses(1,:),'-+black')
plot(losses(2,:),'-+c')
plot(losses(3,:),'-oblack')
plot(losses(4,:),'-oc')
ylim([0 1])
xlim([1 max(k_vec)])
ylabel('RRMSE','FontSize',17);
xlabel('$k$','Interpreter','latex','FontSize',17);
legend(findall(gca,'Tag','Box'), {'symmetric matrix \tau=0.1','random matrix \tau=0.1','symmetric matrix \tau=0.25','random matrix \tau=0.25'});
legend boxoff
box on
