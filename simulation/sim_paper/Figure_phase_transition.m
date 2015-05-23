% Compute recovery probability as function of number of measurements in noiseless settings
% Set simulation parameters
n=150; r=3; noise=0; X_type='low_rank';
tol=10^-5;
RRMSE_thresh=10^-3; % threshold for successful recovery (in relative error)
k_vec = 3:25;
num_sampled_matrices=50;
alg_str = {'svt','optspace','svls'};
measurement_type={'random_entries','random_entries','columns_and_rows'};
alg_iters= [500, 500, 100];

% Output file names
data_file_name = fullfile(master_dir, 'results_vec', ['phase_transition_n_' num2str(n) '_r_' num2str(r) '.m']);
figure_file_name = fullfile(master_dir, 'results_pic', ['phase_transition_n_' num2str(n) '_r_' num2str(r)]);

if(simulate_flag)
    losses_arr = cell(length(alg_str),1); success_arr=losses_arr; times_arr=losses_arr; iters_arr=losses_arr;
    for i=1:length(alg_str)
        [losses_arr{i}, times_arr{i}, iters_arr{i}] = All_Simulations( ...
            n,r,noise,alg_str(i),num_sampled_matrices,alg_iters(i),tol,k_vec ,X_type,measurement_type{i});
        success_arr{i} = losses_arr{i}<RRMSE_thresh;
    end
    save(data_file_name, 'losses_arr', 'success_arr', 'times_arr', 'iters_arr', 'n', 'r', 'k_vec', 'alg_str', ...
        'measurement_type', 'alg_iters', 'num_sampled_matrices', 'tol');
else
    load(data_file_name);
end

% Plot Figure
d_over_n_sqr_vec = 2*n*k_vec/n^2;
symbol_vec = {'s', '*', 'o'};
color_vec = {'g','r','b'};
figure; hold on
for i=1:length(alg_str)
    plot(d_over_n_sqr_vec, mean(success_arr{i}), ['--' color_vec{i} symbol_vec{i}]);
end

ylabel('Reconstruction Rate','FontSize',17)
xlabel('$d/n^2$','Interpreter','latex','FontSize',17);

hLegend = legend(findall(gca,'Tag','Box'), {'SVT MC','optSpace MC','SVLS RCMC'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
set(hChildren(6),'Color','black')
set(hChildren(4),'Color','r')
set(hChildren(2),'Color','b')
legend boxoff

my_saveas(gcf, figure_file_name, {'epsc', 'pdf', 'fig'}); % save figure to file
