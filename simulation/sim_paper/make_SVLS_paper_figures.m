
% Plot all simulation figures shown in the SVLS paper
AssignGeneralConstants;
% Set master directory. Please change this to a directory where SVLS is installed 
if(exist('C:\Users\user\Documents\GitHub\SVLS\', 'dir'))
    master_dir =  'C:\Users\user\Documents\GitHub\SVLS';
else
    master_dir =  'C:\Users\Avishay\Dropbox\matrix_completion\src\my_thesis';
end
figs_dir = fullfile(master_dir, 'figs'); 
data_Dir = fullfile(master_dir, 'data'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters:
simulate_flag = 1; % 1 - run simulations again, 0 - just load results and plot them

figure_type = 'fig';% type of figure files to save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Figure 1: Phase transition plot %%%%%%%%%%%%%%%%%%%%%%
%%%% Probabilty of reconstruction as function of number of measurements in %%%%
%%%% noisless case for different algorithms and designs %%%% %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Figure_phase_transition; % Time Consuming  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2: boxplots for RCMC
%%% Plot the reconstruction error (RMSE) as function of number of
%%% measurements for different algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Figure_Boxplot_RCMC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 3: Compare Gaussian Ensemble and GRC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Figure_APGL_VS_SVLS; % Time Consuming 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Figures 4,5: change_n and change_p_r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Figure_change_n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Figure 6: change tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Figure_change_tau; 
