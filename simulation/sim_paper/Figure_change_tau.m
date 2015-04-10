%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot reconstruction error as function of noise level tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau=[0:0.005:0.1]; % noise level
n=1000; X_type='low_rank';measurement_type= 'columns_and_rows_normal';
k=100; 
num_sampled_matrices=10;
r_vec = [2 4 6 8]; % change rank 

% Output file names 
data_file_name = fullfile(data_dir, ['tau_n_' num2str(n) '_k_' num2str(k) '.m']); 
figure_file_name = fullfile(figs_dir,  ['tau_n_' num2str(n) '_k_' num2str(k) ]); 

if(simulate_flag)
rel_error_mat = zeros(length(tau), length(r_vec)); 
for r = 1:length(r_vec)
for i=1:length(tau)
    for j=1:num_sampled_matrices
        [ X, measure ] = samp_matrix( ...
            n,r_vec(r),k,tau(i), X_type,measurement_type);
        X_hat = SVLS( measure.Br, measure.Bc, measure.Ar, measure.Ac,[]);
        rel_error_mat(i,r) = rel_error_mat(i,r)  + RRMSE(X_hat, X); % add error (we will divide by T  in the end) %  +norm(m-N,'fro')/(T*norm(m,'fro'));
    end
end
end
rel_error_mat = rel_error_mat ./ num_sampled_matrices; % normalize by number of iterations 
	    save(data_file_name, 'rel_error_mat', 'num_sampled_matrices', 'n', 'k', 'r_vec', ...
        'measurement_type', 'num_sampled_matrices');
else
	load(data_file_name); 
end

%%%Plot Figure
figure; hold on
color_vec={'r','g','k','m'};
%color_vec = {'r','b','g','m'}
for r=1:length(r_vec)
	plot(rel_error_mat(:,r), color_vec{r}); 
end

ylim([0 0.04]);
xlabel('$\tau$','Interpreter','latex','FontSize',17);   ylabel('RRMSE','FontSize',17);
hLegend = legend(findall(gca,'Tag','Box'), {'r=2','r=4','r=6','r=8'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
legend boxoff

my_saveas(gcf, figure_file_name, {'epsc', 'pdf', 'fig'}); % save figure to file 
