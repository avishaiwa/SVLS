% Figure 2: boxplots for RCMC
%%% Plot the reconstruction error (RMSE) as function of number of
%%% measurements for different algorithms
n=150; r=3;
noise=0.25; % This is sigma or sigma^2?
X_type=  'low_rank';measurement_type= 'columns_and_rows';
alg_str = {'svls','optspacenotrim','svt_col_row'};
num_sampled_matrices=50;
alg_iters=100; tol=10^-5;
k_vec= (8:2:16);

% Output file names
DataFileName = fullfile(data_dir, ['boxplot_n_' num2str(n) '_r_' num2str(r) '_sigma_' num2str(noise) '.mat']);
FigureFileName = fullfile(figs_dir, ['boxplot_n_' num2str(n) '_r_' num2str(r) '_sigma_' num2str(noise) '.mat']);
%data_file_name = fullfile(master_dir, 'results_vec', ['boxplot_n_' num2str(n) '_r_' num2str(r) '_sigma_' num2str(noise) '.mat']);
%figure_file_name = fullfile(master_dir, 'results_pic', ['boxplot_n_' num2str(n) '_r_' num2str(r)  '_sigma_' num2str(noise) ]);

if(simulate_flag)
    [ losses, times, iterations ] = All_Simulations( ...
        n, r, noise, alg_str, num_sampled_matrices, alg_iters, tol, k_vec ,X_type, measurement_type);
    save(fullFileName, 'losses', 'times', 'iterations', 'n', 'r', 'k_vec', 'alg_str', ...
        'measurement_type', 'alg_iters', 'num_sampled_matrices', 'tol');
else
    load(data_file_name);
end

% Make figure
figure; hold on
for i=1:length(alg_str)
    if(i==2)
        x_lab = num2str_cell(k_vec');
    else
        x_lab = repmat({''}, 1, length(k_vec));
    end
    color_vec = ['b','r','g'];
    boxplot(log(losses(:,:,i)), 'color', color_vec(i), 'symbol','+', ...
        'positions', 0.25*(i-2) + (1:5), 'labels', x_lab, 'widths',0.2);
end

ylim([min(log(losses(:))) max(log(losses(:)))]);  % ylim([-2.5 1.5]);
ylabel('RRMSE','FontSize',17);
xlabel('$k$','Interpreter','latex','FontSize',17);

hLegend = legend(findall(gca,'Tag','Box'), {'SVT','optSpace','SVLS'}); % why is this reversed????
hChildren = findall(get(hLegend,'Children'), 'Type','Line'); % Why is this needed?
set(hChildren(6),'Color','g');
set(hChildren(4),'Color','r');
set(hChildren(2),'Color','b');
legend boxoff;

my_saveas(gcf, FigureFileName, {'epsc', 'pdf', 'fig'}); % save figure to file
