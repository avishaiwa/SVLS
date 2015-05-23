if(exist('C:\Users\user\Dropbox\matrix completion\', 'dir'))
    master_dir =  'C:\Users\user\Dropbox\matrix completion\src\my_thesis';
else
     master_dir =  'C:\Users\Avishay\Dropbox\matrix completion\src\my_thesis';
end

% output_fig_name = [sampling_scheme '_' alg_str];
% output_figs_dir = 'figs/';
% my_saveas(gcf, fullfile(output_figs_dir, output_fig_name), {'fig', 'epsc'}); % save in different formats
% 
% my_saveas(gcf, fullfile(master_dir, 'results_pic', '100_2_NEW'), {'epsc', 'pdf', 'jpg'});
% 
%%%%%
%%%%% phase transition plot
figure;
load(fullfile(master_dir, 'results_vec2', 'phase_transition_150_3.mat'));
plot(index,result_SVT,'--sk')
hold on
plot(index,result_opt,'--r*')
plot(index,result_SVLS,'--o')

ylabel('Reconstruction Rate','FontSize',17)
xlabel('$d/n^2$','Interpreter','latex','FontSize',17);

hLegend = legend(findall(gca,'Tag','Box'), {'SVT MC','optSpace MC','SVLS RCMC'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
set(hChildren(4),'Color',[0 0 1])
set(hChildren(6),'Color','g')
set(hChildren(4),'Color','r')
set(hChildren(2),'Color','b')
legend boxoff

%%%%%%

% boxplots for RCMC 
figure;
load(fullfile(master_dir, 'results_vec2', 'RCMC_100_5_025.mat'));

boxplot(log(losses(:,:,1)),'COLOR','r','symbol','+','positions',[0.75 1.75 2.75 3.75 4.75],'labels',{'','','','',''},'widths',0.2)
hold on
boxplot(log(losses(:,:,2)),'COLOR','k','symbol','+','positions',[1.25 2.25 3.25 4.25 5.25],'labels',{'','','','',''},'widths',0.2)
boxplot(log(losses(:,:,3)),'COLOR','b','symbol','+','labels',{'8','10','12','14','16'},'widths',0.2)
ylim([-2.5 1.5])

ylabel('RRMSE','FontSize',17)
xlabel('$k$','Interpreter','latex','FontSize',17);

hLegend = legend(findall(gca,'Tag','Box'), {'SVT','optSpace','SVLS'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
set(hChildren(6),'Color','k')
set(hChildren(4),'Color','r')
set(hChildren(2),'Color','b')
legend boxoff

%%%%%%%%%%%%%%%%%%%%
%%%gaussian and GRC 
figure;
load(fullfile(master_dir, 'results_vec2', 'APGL_basis_100_2.mat'));
loglog_flag=1;
if(loglog_flag)
    loglog(index,lossesG01,'-.^b')
    
    hold on
    loglog(index,lossesG001,'--+b')
    loglog(index,lossesG0001,'-*b')
    
    loglog(index,lossesCR01,'-.sr')
    loglog(index,lossesCR001,'--dr')
    loglog(index,lossesCR0001,'-or')
else
    semilogy(log(index),lossesG01,'-.^b')
    
    hold on
    semilogy(log(index),lossesG001,'--+b')
    semilogy(log(index),lossesG0001,'-*b')
    
    semilogy(log(index),lossesCR01,'-.sr')
    semilogy(log(index),lossesCR001,'--dr')
    semilogy(log(index),lossesCR0001,'-or')
end

xlim([0.038 0.4]); % xlim([10^-1.5 3])
set(gca,'XTick',(index))
ylabel('RRMSE','FontSize',17)
xlabel('$d/n^2$','Interpreter','latex','FontSize',17);
hLegend = legend(findall(gca,'Tag','Box'), {'GE \tau=0.1','GE \tau=0.01','GE \tau=0.001','GRC \tau=0.1','GRC \tau=0.01','GRC \tau=0.001'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');


legend boxoff

%%%%%%%%%%%%%%%%%%%%
%%%changen_p_r

figure;
load(fullfile(master_dir, 'results_vec2', 'change_n_p_r_.mat'));
index=[20:20:20*50];

plot(index,mean(result'))
hold on
plot(index,mean(result0'),'-r')
ylabel('RRMSE','FontSize',17)
xlabel('$n$','Interpreter','latex','FontSize',17);
hLegend = legend(findall(gca,'Tag','Box'), {' \tau=0.1','\tau=0.01'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
legend boxoff

%%%%
%%%change_n
figure;
load(fullfile(master_dir, 'results_vec2', 'change_n.mat'));
index=[10:10:1000];
plot(index,mean(result1'))
hold on
plot(index,mean(result'),'-r')
ylabel('RRMSE','FontSize',17)
xlabel('$n$','Interpreter','latex','FontSize',17);
hLegend = legend(findall(gca,'Tag','Box'), {' \tau=0.1','\tau=0.01'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
legend boxoff

%%%%
%%%change tau
figure;
load(fullfile(master_dir, 'results_vec2', 'tau_1000_100.mat'));
plot(result2,'-r');
hold on
plot(result4,'-b'); plot(result6,'-g'); plot(result8,'-k'); ylim([0 0.04])
xlabel('$\tau$','Interpreter','latex','FontSize',17);   ylabel('RRMSE','FontSize',17);
hLegend = legend(findall(gca,'Tag','Box'), {'r=2','r=4','r=6','r=8'});
hChildren = findall(get(hLegend,'Children'), 'Type','Line');
legend boxoff


