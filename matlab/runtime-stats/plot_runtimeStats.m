%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Runtime Computation plots for SVM models of increasing
% complexity, models: 
% 0: 1.6kSVs (0.5% Data),1: 2.7kSVs (1% Data), 2: 3.5kSVs (2% Data)
% 3: 8.9kSVs (5%  Data), 4: 15kSVs  (10% Data) from 1 million points.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

%% Extract Data from Computation time files
model_folders = getfolders('.');

mean_coll = zeros(1,length(model_folders));
std_coll  = zeros(1,length(model_folders));

mean_IK = zeros(1,length(model_folders));
std_IK  = zeros(1,length(model_folders));

for i=1:length(model_folders)   
    clear coll_times IK_times
    coll_times = textread(strcat('./',model_folders(i).name,'/Collision_boundary_construction_time.txt'), '%f')*1000;
    IK_times   = textread(strcat('./',model_folders(i).name,'/IK_solver_computation_time.txt'), '%f')*1000;

    % Better than mean() and std()
    [muhat_coll,sigmahat_coll] = normfit(coll_times);
    [muhat_IK,  sigmahat_IK]   = normfit(IK_times);
    
    mean_coll(i) = muhat_coll; std_coll(i) = sigmahat_coll;
    mean_IK(i)   = muhat_IK;   std_IK(i)   = sigmahat_IK;
    
end

model_size  =  [1600 2700 3500 8900 15000];
data_perc   =  [0.5  1     2    5 10 ];

%% Plot Computation Times
clc; close all;
figure('Color', [1 1 1]);
errorbar(model_size,mean_coll, std_coll,'--or','LineWidth',2); hold on;
errorbar(model_size,mean_IK, std_IK,'--ob','LineWidth',2); hold on;
% plot(model_size,mean_IK,'--ob','LineWidth',2); hold on;
set(gca,'xscale','log', 'XMinorTick','on', 'XMinorGrid','on')
xlim([1000 20000])
grid on
xlabel('$N_{sv}$','Interpreter','Latex','FontSize',14); ylabel('Computation time [ms]','Interpreter','Latex','FontSize',14)
legend({'$\Gamma(x)$ and $\nabla\Gamma(x)$ Evaluation ', 'IK Solver'},'Interpreter','Latex', 'FontSize',14)

%% Plot Classifier Performance
figure('Color',[1 1 1])
[ax,hline1,hline2]=plotyy(model_size, mean_coll,model_size,data_perc);
delete(hline1);
delete(hline2);
hold on;
plot(model_size,data_perc, '--og','LineWidth',2)
set(gca,'xscale','log')
xlabel('$N_{sv}$','Interpreter','Latex','FontSize',14);
