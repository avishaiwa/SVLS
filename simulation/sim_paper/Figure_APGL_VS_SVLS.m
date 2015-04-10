%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 3: Compare Gaussian Ensemble and GRC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set parameters
n=100; r=2; X_type='low_rank'; 
k_vec=[2:10,12:2:20]; % number of row and columns measurements 
num_sampled_matrices = 10; % how many samples for each parameter value

tau_vec = [0.001 0.01 0.1]; % vector of noise levels
alg_str = {'SVLS', 'apgl'}; % vector of algorithms
losses_array = cell(length(tau_vec), length(alg_str)); 
times_array = losses_array; iters_array = losses_array;

% Output file names 
DataFileName = fullfile(data_dir, ['boxplot_n_' num2str(n) '_r_' num2str(r) '_sigma_' num2str(noise) '.mat']);
FigureFileName = fullfile(figs_dir, ['boxplot_n_' num2str(n) '_r_' num2str(r) '_sigma_' num2str(noise) '.mat']);
%data_file_name = fullfile(master_dir, 'results_vec', ['boxplot_n_' num2str(n) '_r_' num2str(r) '_sigma_' num2str(noise) '.mat']);
%figure_file_name = fullfile(master_dir, 'results_pic',%['boxplot_n_'%num2str(n) '_r_' num2str(r)  '_sigma_' num2str(noise) ]);

for i=1:length(tau_vec) % run and compute losses
    for j=1:length(alg_str)
        if(strcmp(alg_str{j}, 'SVLS'))
            alg_num_iters=100; % set # iterations
            tol = 10^-5; 
            measurement_type= 'columns_and_rows_normal';
        else % apgl
            alg_num_iters=2000; 
            tol = 10^-7;
            measurement_type= 'gaussian';
        end
        [losses_array{i,j}, times_array{i,j}, iters_array{i,j} ] = All_Simulations( ...
            n,r,tau_vec(i),alg_str(j),num_sampled_matrices,alg_num_iters,tol,k_vec ,X_type,measurement_type);
    end
end
% save parameters and results 
save(DataFileName, 'losses_array', 'times_array', 'iters_array', ...
    'n', 'r', 'X_type', 'k_vec', 'tau_vec', 'alg_str'); 


% .... (missing-here) % Save results to file 


% % % %%SVLS with tau=0.001
% % % n=100;r=2;noise=0.001;X_type=  'low_rank';measurement_type= 'columns_and_rows_normal';
% % % alg_str = {'SVLS'}; p_iter=50; num_iter=100; tol=10^-5; index=[2:10,12:2:20];
% % % 
% % % [ lossesCR0001 times, iterations ] = All_Simulations( ...
% % %     n,r,noise,alg_str,p_iter,num_iter,tol,index ,X_type,measurement_type);
% % % 
% % % %%APGL with tau=0.001
% % % n=100;r=2;noise=0.001;X_type=  'low_rank';measurement_type= 'gaussian';
% % % alg_str = {'apgl'}; p_iter=10; num_iter=2000; tol=10^-7; index=[2:10,12:2:20];
% % % 
% % % [ lossesG0001, times, iterations ] = All_Simulations( ...
% % %     n,r,noise,alg_str,p_iter,num_iter,tol,index ,X_type,measurement_type);
% % % 
% % % %%SVLS with tau=0.01
% % % n=100;r=2;noise=0.01;X_type=  'low_rank';measurement_type= 'columns_and_rows_normal';
% % % alg_str = {'SVLS'}; p_iter=50; num_iter=100; tol=10^-5; index=[2:10,12:2:20];
% % % 
% % % [ lossesCR001 times, iterations ] = All_Simulations( ...
% % %     n,r,noise,alg_str,p_iter,num_iter,tol,index ,X_type,measurement_type);
% % % 
% % % %%APGL with tau=0.01
% % % n=100;r=2;noise=0.01;X_type=  'low_rank';measurement_type= 'gaussian';
% % % alg_str = {'apgl'}; p_iter=10; num_iter=2000; tol=10^-7; index=[2:10,12:2:20];
% % % 
% % % [ lossesG001, times, iterations ] = All_Simulations( ...
% % %     n,r,noise,alg_str,p_iter,num_iter,tol,index ,X_type,measurement_type);
% % % 
% % % %%SVLS with tau=0.1
% % % n=100;r=2;noise=0.1;X_type=  'low_rank';measurement_type= 'columns_and_rows_normal';
% % % alg_str = {'SVLS'}; p_iter=50; num_iter=100; tol=10^-5; index=[2:10,12:2:20];
% % % 
% % % [ lossesCR01 times, iterations ] = All_Simulations( ...
% % %     n,r,noise,alg_str,p_iter,num_iter,tol,index ,X_type,measurement_type);
% % % 
% % % %%APGL with tau=0.1
% % % n=100;r=2;noise=0.1;X_type=  'low_rank';measurement_type= 'gaussian';
% % % alg_str = {'apgl'}; p_iter=10; num_iter=2000; tol=10^-7; index=[2:10,12:2:20];
% % % 
% % % [ lossesG01, times, iterations ] = All_Simulations( ...
% % %     n,r,noise,alg_str,p_iter,num_iter,tol,index ,X_type,measurement_type);


%Make Figure
d_over_n_square_vec=2*n*k_vec/n^2; % transfer k to d/n^2
figure;      hold on

plot_symbols_vec = {'-*b', '--+r', '-.^b', '-or', '--db', '-.sr'}; 
for i=1:length(tau_vec) % run and compute losses
    for j=1:length(alg_str)
        loglog(d_over_n_square_vec, mean(losses_array{i,j}), [ plot_symbols_vec{2*(i-1)+j} ]);
% loglog(index, mean(lossesG001),'--+r');
% loglog(index, mean(lossesG0001), '-*r');
% loglog(index, mean(lossesCR01), '-.sb');
% loglog(index, mean(lossesCR001), '--db');
% loglog(index, mean(lossesCR0001), '-ob');
    end
end
xlim([0.038 0.4]); % xlim([10^-1.5 3])
set(gca,'XTick', d_over_n_square_vec);
ylabel('RRMSE','FontSize',17);
xlabel('$d/n^2$','Interpreter','latex','FontSize',17);
hLegend = legend(findall(gca,'Tag','Box'), ...
    {'GRC \tau=0.001','GE \tau=0.001','GRC \tau=0.01','GE \tau=0.01','GRC \tau=0.1','GE \tau=0.001'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
legend boxoff
my_saveas(gcf, FigureFileName, {'epsc', 'pdf', 'fig'}); % save figure to file 
