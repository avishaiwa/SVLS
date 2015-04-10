%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Two Figures: change_n k and r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First Figure: tau=0.01, changing n
r=3;
tau_vec=[0.01 0.1];
X_type='low_rank';measurement_type= 'columns_and_rows_normal';
k=20;
num_sampled_matrices=10; % number of iterations
n_vec = 10:10:1000; % set n
color_vec={'b','r'};


% Output file names
DataFileName = fullfile(data_dir, ['boxplot_n_' num2str(n) '_r_' num2str(r) '_sigma_' num2str(noise) '.mat']);
FigureFileName = fullfile(figs_dir, ['boxplot_n_' num2str(n) '_r_' num2str(r) '_sigma_' num2str(noise) '.mat']);
%data_file_name = fullfile(master_dir, 'results_vec', ['boxplot_n_' num2str(n) '_r_' num2str(r) '_sigma_' num2str(noise) '.mat']);
%figure_file_name = fullfile(master_dir, 'results_pic', ['boxplot_n_'%num2str(n) '_r_' num2str(r)  '_sigma_' num2str(noise) ]);

if(simulate_flag)
    rel_error_mat = zeros(length(n_vec),num_sampled_matrices,length(tau_vec));
    for t=1:length(tau_vec)
        for i=1:length(n_vec)
            for j=1:num_sampled_matrices
                [ X, measure ] = samp_matrix( ...
                    n_vec(i),r,k,tau_vec(t), X_type,measurement_type);
                X_hat = SVLS( measure.Br,measure.Bc,measure.Ar,measure.Ac,[]);
                rel_error_mat(i,j,t)= RRMSE(X_hat, X); %  norm(m-N,'fro')/(norm(m,'fro'))
            end
        end
    end
    
    save(DataFileName , 'rel_error_mat', 'n_vec', 'r', 'k', 'tau_vec', ...
        'measurement_type', 'num_sampled_matrices', 'tol');
else
    load(data_file_name);
end

%% Plot figure
figure; hold on
for t=1:length(tau_vec)
    plot(n_vec, mean(rel_error_mat(:,:,t),2), color_vec{t});
end
ylabel('RRMSE','FontSize',17)
xlabel('$n$','Interpreter','latex','FontSize',17);
hLegend = legend(findall(gca,'Tag','Box'), {' \tau=0.1','\tau=0.01'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
legend boxoff
my_saveas(gcf, FigureFileName,  {'epsc', 'pdf', 'fig'});%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Second Figure: tau=0.01, changing n,k and r

if(simulate_flag)
    rel_error_mat = zeros(length(n_vec),num_sampled_matrices,length(tau_vec));
    r_vec = 1:100; n_vec = 20.*r_vec; k_vec=4.*r_vec; % set n,k,r
    data_file_name = fullfile(master_dir, 'results_vec2', ['change_n_p_r_oversamp_4_n_over_r_20' '.m']);
    figure_file_name = fullfile(master_dir, 'results_pic2', ['change_n_p_r_oversamp_4_n_over_r_20' ]);
    
    
    for t=1:length(tau_vec)
        for i=1:length(n_vec)
            for j=1:num_sampled_matrices
                [ X, measure ] = samp_matrix( ...
                    n_vec(i),r_vec(i),k_vec(i),tau_vec(t), X_type,measurement_type  );
                X_hat = SVLS( measure.Br,measure.Bc,measure.Ar,measure.Ac,[]);
                rel_error_mat(i,j,t) = RRMSE(X_hat,X);
            end
        end
    end
    save(DataFileName , 'rel_error_mat', 'n_vec', 'r_vec', 'k_vec', 'tau_vec', ...
        'measurement_type', 'num_sampled_matrices', 'tol');
else
    load(data_file_name);
end

%% Plot figure
figure; hold on
for t=1:length(tau_vec)
    plot(n_vec, mean(rel_error_mat(:,:,t),2), color_vec{t});
end

ylabel('RRMSE','FontSize',17);
xlabel('$n$','Interpreter','latex','FontSize',17);
hLegend = legend(findall(gca,'Tag','Box'), {' \tau=0.1','\tau=0.01'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
legend boxoff
my_saveas(gcf, FigureFileName,  {'epsc', 'pdf', 'fig'});%%%%%%
