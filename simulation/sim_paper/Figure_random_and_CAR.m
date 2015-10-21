n=500; r=10; noise=0.5;alg_str={'SVLS_sym'};num_sampled_matrices=100; measurement_type='columns_and_rows'; k =50;
max_iter=50;tol=10^-4;    X_type='low_rank';
num=5;
result1 = zeros(1,21);
result_opt = zeros(1,21);
j=0;
for p=0:1000:20000
    p
    j=j+1;
    for i=1:num
        [ m, measure ] = samp_matrix( ...
            n,r,k,noise, X_type,measurement_type  );
        Br =  measure.Br; Bc =  measure.Bc; Ar =  measure.Ar; Ac = measure.Ac;
        [n1 n2 ] = size(Bc);
        aBc = Bc;
        samp=randsample(n1*n2,p); aBc(samp) = 0;
        [X S Y dist] = OptSpace(aBc,r,max_iter,tol);
        aBc1 = X*S*Y';
        aBc1(aBc~=0) = Bc(aBc~=0);
        
        [n1 n2 ] = size(Br);
        aBr = Br;
        samp=randsample(n1*n2,p); aBr(samp) = 0;
        [X S Y dist] = OptSpace(aBr,r,max_iter,tol);
        aBr1 = X*S*Y';
        aBr1(aBr==0) = Br(aBr==0);
        
        [X_hat] = SVLS_sym( aBr1, aBc1, Ar, Ac,r);
        RRMSE(X_hat,m)
        result1(j) = result1(j)  +RRMSE(X_hat,m)/num;
        
        am = measure.am;
        samp = randsample(find(am~=0),2*p);
        for s=samp
            am(s)=0;
        end
        [X S Y dist] =OptSpacenoTrim(am,r,max_iter,tol);
        X_hat =  X*S*Y';
        RRMSE(X_hat,m)
        result_opt(j) = result_opt(j)  +RRMSE(X_hat,m)/num;
        
    end
end

plot(1-(0:1000:20000)/25000,result1)
hold on
plot(1-(0:1000:20000)/25000,result_opt,'-r')
xlabel('p','FontSize',17)
ylabel('RRMSE','FontSize',17)
legend(findall(gca,'Tag','Box'), {'SVLS+optSpace','optSpace'});
legend boxoff
box on
