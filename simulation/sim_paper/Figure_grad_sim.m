n=[200]; r=4; noise=0.5;alg_str={'SVLS_sym'};num_sampled_matrices=100; measurement_type='columns_and_rows'; k =40;
max_iter=30;tol=10^-100;    X_type='low_rank';


[ m, measure ] = samp_matrix( ...
                n,r,k,noise, X_type,measurement_type  );
            
[X_hat] = SVLS_sym( measure.Br, measure.Bc, measure.Ar, measure.Ac,[]);

tic
[X,dis,result] = gradientDescent( X_hat, measure.Br, measure.Bc, measure.Ar, measure.Ac,r, ...
    'gradient_descent'  ,tol, max_iter,m );
t1 = toc;
tic
[X_rand,dis_rand,result_rand] = gradientDescent( randn(n), measure.Br, measure.Bc, measure.Ar, measure.Ac,r, ...
    'gradient_descent'  ,tol, max_iter,m );
t2 = toc
tic
 [X1,dis1,result1] = Alternating_LS(X_hat, measure.Br, measure.Bc, measure.Ar, measure.Ac,r, ...
    tol, max_iter,1,m );
t3 = toc;
tic
 [X1_rand,dis1_rand,result1_rand] = Alternating_LS3(randn(n), measure.Br, measure.Bc, measure.Ar, measure.Ac,r, ...
    tol, max_iter,1,m );
t4 = toc;
%result_rand(result_rand>1)=1;
%result1_rand(result1_rand>1)=1;
%result1(result1>1)=1;


figure
plot(result)
xlabel('number of iterations')
ylabel('RRMSE')
hold on
plot(result1,'-r'); plot(result_rand,'-black'); plot(result1_rand,'-g'); 
legend(findall(gca,'Tag','Box'), {'gradient descent','LS'});
legend boxoff
box on
figure
loglog(dis)
xlabel('number of iterations')
ylabel('||A^{(R)}X-B^{(R)}||_F+||XA^{(C)}-B^{(C)}||_F')
hold on; loglog(dis1,'col','red'); plot(dis_rand,'-black'); plot(dis1_rand,'-g'); 
legend(findall(gca,'Tag','Box'), {'gradient descent','LS'});
legend boxoff
box on

figure
plot(result)
xlabel('number of iterations','FontSize',17)
ylabel('RRMSE','FontSize',17)
hold on
plot(result1,'-r');
legend(findall(gca,'Tag','Box'), {'gradient descent','LS'});
legend boxoff
box on

figure
loglog(dis)
xlabel('number of iterations','FontSize',17)
ylabel('||A^{(R)}X-B^{(R)}||_F+||XA^{(C)}-B^{(C)}||_F','FontSize',17)
hold on; loglog(dis1,'col','red');
legend(findall(gca,'Tag','Box'), {'gradient descent','LS'});
legend boxoff
box on