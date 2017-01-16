%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training a model with SVM-Light and creating an SVMGrad Data-structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load 36D Robot Self-Collision Dataset and model learned through libSVM
clc; clear all; close all;
load('./models/36d-robot-collision.mat')

%% Create Simplified Struct Model for SVMGrad from libSVM Model
svmgrad = [];
svmgrad.D       = size(X_train,2);
svmgrad.nSV     = model.totalSV;
svmgrad.b       = -model.rho;
svmgrad.sigma   = options.sigma;
svmgrad.yalphas = model.sv_coef'; %\alpha_*y_i
svmgrad.SVs     = full(model.SVs)';

%% Example of learning an C-SVM with SVMLight

% Input Parameters
% -t : kernel_type (2:RBF)
% -g : gamma (gamma = 1/2sigma^2))
% -C : Misclassifiction penalty

C     = 223;
sigma = 0.633;
gamma = 1/(2*svmgrad.sigma*svmgrad.sigma);
tic;
model = svmlearn(X_train, y_train,'-v 3 -z c -c 223 -t 2 -g 1.2479 ');
toc;