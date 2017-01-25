%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Performance Computation plots for SVM models of increasing
% complexity, models: 
% 0: 1.6kSVs (0.5% Data),1: 2.7kSVs (1% Data), 2: 3.5kSVs (2% Data)
% 3: 8.9kSVs (5%  Data), 4: 15kSVs  (10% Data) from 1 million points.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;
%% Load full dataset (Train and Test)
load('./models/Full-Collision-Avoidance-Dataset.mat')

%% Load models
model_names = dir('./models/*36d-*.mat'); 
nModels = length(model_names);
models  = {};
for i=1:nModels
    models{i,1} = load(strcat('./models/',model_names(i).name));
end
%% Compute Error Rates for each model on test set
numSamples  = round(length(y_test));
model_sizes = zeros(1,nModels);    
ACC = zeros(1,nModels);
F1  = zeros(1,nModels);
FPR = zeros(1,nModels);
TPR = zeros(1,nModels);
clc;
for i =  1:nModels    
    test_model = models{i}.model;
    model_sizes(1,i) = models{i}.model.totalSV;
    
    % Extract a random testing points
    y_test_   = zeros(1,numSamples);
    y_est_    = zeros(1,numSamples);
    idx_rand  = randperm(numSamples);
    
    % Evaluate on Testing Points
    for ii=1:numSamples
        X_test_    = X_test(idx_rand(ii),:);
        y_test_(ii) = y_test(idx_rand(ii),:);
        
        % Test Learnt Model
        [y_est_(ii)] = svm_classifier(X_test_, y_test_(ii), [], test_model);
    end
    
    % Compute Classifier Test Stats
    [stats] = class_performance(y_test_,y_est_);
    fprintf('*Classifier Performance for N_sv=%d on Test Set (%d points)* \n Acc: %1.5f, F-1: %1.5f, FPR: %1.5f, TPR: %1.5f \n', test_model.totalSV, length(y_test_), stats.ACC, stats.F1, stats.FPR, stats.TPR)
    ACC(1,i) = stats.ACC; F1(1,i) = stats.F1;
    FPR(1,i) = stats.FPR; TPR(1,i) = stats.TPR;
end

%% Plot of model performance stats on test set
figure('Color',[1 1 1])
plot(model_sizes, ACC, '--og','LineWidth',2); hold on;
plot(model_sizes, F1,  '--vk','LineWidth',2); hold on;
plot(model_sizes, 1- FPR, '--dr','LineWidth',2); hold on;
plot(model_sizes, TPR, '--sb','LineWidth',2); 
xlabel('$N_{sv}$','Interpreter','Latex','FontSize',14);
ylabel('Performance Measure','Interpreter','Latex','FontSize',14)
set(gca,'xscale','log')
xlim([1000 20000])
legend({'ACC','F1','1-FPR','TPR'},'Interpreter','Latex', 'FontSize',14)
grid on
