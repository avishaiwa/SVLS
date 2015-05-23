%%rank estimation
n=400; r=7 ;alg_str={'svls'};num_sampled_matrices=10; measurement_type='columns_and_rows'; k = 15;X_type='low_rank';
maximum_iter=10;tol=10^-4;long=100;noise=0.25;


% Output file names


k_vec7 =  (1:1:60);
rank_est7=zeros(1,length(k_vec7));


for i=2:length(k_vec7)
    for j=1:num_sampled_matrices
        [ m, measure ] = samp_matrix( ...
            n,r,k_vec7(i),noise, X_type,measurement_type  );
        [X_hat, r_est] = SVLS( measure.Br, measure.Bc, measure.Ar, measure.Ac,[]);
        rank_est7(i)=rank_est7(i)+r_est/num_sampled_matrices;
    end
end

r=15
k_vec15 = (1:1:60);
rank_est15=zeros(1,length(k_vec15));
rank_est15(1:30)=0;
for i=2:length(k_vec15)
    for j=1:num_sampled_matrices
        [ m, measure ] = samp_matrix( ...
            n,r,k_vec15(i),noise, X_type,measurement_type  );
        [X_hat, r_est] = SVLS( measure.Br, measure.Bc, measure.Ar, measure.Ac,[]);
        rank_est15(i)=rank_est15(i)+r_est/num_sampled_matrices;
    end
end


num_sampled_matrices=100;
r=7
k_vec7 =  (1:1:60);
rank_est_opt7=7*ones(1,length(k_vec7));


for i=56:60
    i
    rank_est_opt7(i) = 0
    for j=1:num_sampled_matrices
        [ m, measure ] = samp_matrix( ...
            n,r,k_vec7(i),noise, X_type,'random_entries'  );
 [X S Y dist] = OptSpace(measure.am,[],1,1);
        rank_est_opt7(i)=rank_est_opt7(i)+rank(S)/num_sampled_matrices;
    end
end

r=15
k_vec15 =  (1:1:60);
rank_est_opt15=zeros(1,length(k_vec7));

for i=1:60
    i
    for j=1:num_sampled_matrices
        [ m, measure ] = samp_matrix( ...
            n,r,k_vec7(i),noise, X_type,'random_entries'  );
 [X S Y dist] = OptSpace(measure.am,[],1,1);
        rank_est_opt15(i)=rank_est_opt15(i)+rank(S)/num_sampled_matrices;
    end
end
 data_file_name = fullfile(['results_vec\est_rank' '.mat']);
 figure_file_name = fullfile(master_dir, 'results_pic\rank_est');

 
save( data_file_name  , 'rank_est7', 'rank_est_opt7', 'rank_est15','rank_est_opt15');
 d_vec= 2*n*k_vec15-k_vec15.^2;
plot( d_vec,rank_est7,'-b'); hold on; plot( d_vec,rank_est15,'--b');  plot( d_vec,rank_est_opt7,'-r');  plot(d_vec,rank_est_opt15,'--r'); 
xlabel('$d$','Interpreter','latex','FontSize',17);
ylabel('Estimation Rank','FontSize',17);
set(gca, 'ylim', [0 17.5]);
hLegend = legend(findall(gca,'Tag','Box'), {' RCMC, r=7','RCMC, r=15', 'MC, r=7','MC, r=15' });
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
legend boxoff

my_saveas(gcf, figure_file_name, {'epsc', 'pdf', 'fig'}); % save figure to file

