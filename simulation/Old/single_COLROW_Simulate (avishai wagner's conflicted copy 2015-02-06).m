        
        % Collect and Plot results. Need to also save data to file (to get tables)
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

plot(index,result_SVT,'--sk')
hold on
plot(index,result_opt,'--r*')
plot(index,result_basis,'--o')
        ylabel('Reconstruction Rate','FontSize',17)
        xlabel('$d/n^2$','Interpreter','latex','FontSize',17);

        hLegend = legend(findall(gca,'Tag','Box'), {'SVT MC','optSpace MC','SVLS RCMC'});
        hChildren = findall(get(hLegend,'Children'), 'Type','Line');
        set(hChildren(4),'Color',[0 0 1])
         set(hChildren(6),'Color','g')
       set(hChildren(4),'Color','r')
       set(hChildren(2),'Color','b')
       legend boxoff
        
        output_fig_name = [sampling_scheme '_' alg_str];
        output_figs_dir = 'figs/';
        my_saveas(gcf, fullfile(output_figs_dir, output_fig_name), {'fig', 'epsc'}); % save in different formats
        
                      semilogy((2*n*index),(lossesG01),'-.^b')

     hold on
              semilogy((2*n*index),((lossesG001)),'--+b')
                            semilogy((2*n*index),((lossesG0001)),'-*b')

              semilogy((2*n*index),(mean((lossesCR01))),'-.sr')
             semilogy((2*n*index),(mean(lossesCR001)),'--dr')
              semilogy((2*n*index),(mean(lossesCR0001)),'-or')

  ylabel('rel.error','FontSize',17)
        xlabel('$d$','Interpreter','latex','FontSize',17);
        hLegend = legend(findall(gca,'Tag','Box'), {'GE \tau=0.1','GE \tau=0.01','GE \tau=0.001','GRC \tau=0.1','GRC \tau=0.01','GRC \tau=0.001'});
        hChildren = findall(get(hLegend,'Children'), 'Type','Line');
       
       legend boxoff
       
       
       loglog((2*n*index)/n^2,(lossesG01),'-.^b')
       
       hold on
       loglog((2*n*index)/n^2,(lossesG001),'--+b')
       loglog((2*n*index)/n^2,((lossesG0001)),'-*b')
       
       loglog((2*n*index)/n^2,mean((lossesCR01)),'-.sr')
       loglog((2*n*index)/n^2,mean((lossesCR001)),'--dr')
       loglog((2*n*index)/n^2,mean((lossesCR0001)),'-or')
       ylabel('RRMSE','FontSize',17)
       xlabel('$d/n^2$','Interpreter','latex','FontSize',17);
       hLegend = legend(findall(gca,'Tag','Box'), {'GE \tau=0.1','GE \tau=0.01','GE \tau=0.001','GRC \tau=0.1','GRC \tau=0.01','GRC \tau=0.001'});
       hChildren = findall(get(hLegend,'Children'), 'Type','Line');
       
       legend boxoff
       set(gca,'XTick',[0.04:0.04:0.12,0.2,0.3])
       xlim([0.038 0.4])
       saveas(gcf, 'APGL_100_2', 'epsc'); 
       

       
       save('APGL_basis_100_2_two','lossesG01','lossesG001','lossesG0001','lossesCR01','lossesCR001','lossesCR0001','index')
       
       
plot(index,mean(result1'))
hold on
plot(index,mean(result'),'-r')
  ylabel('RRMSE','FontSize',17)
        xlabel('$n$','Interpreter','latex','FontSize',17);
        hLegend = legend(findall(gca,'Tag','Box'), {' \tau=0.1','\tau=0.01'});
        hChildren = findall(get(hLegend,'Children'), 'Type','Line');
       index=[20:20:20*50]
       legend boxoff
        save('change_n','result','result1','index')
        
        
        plot(result2,'-r'); 
        hold on
        plot(result4,'-b'); plot(result6,'-g'); plot(result8,'-k'); ylim([0 0.04])
xlabel('$\tau$','Interpreter','latex','FontSize',17);   ylabel('RRMSE','FontSize',17);
        hLegend = legend(findall(gca,'Tag','Box'), {'r=2','r=4','r=6','r=8'});
        hChildren = findall(get(hLegend,'Children'), 'Type','Line');
              legend boxoff

    
        