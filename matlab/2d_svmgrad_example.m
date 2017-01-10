%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SVMgrad is a compact library used to evaluate the decision function of a
% Gaussian RBF Kernel Support Vector Machine, as well as the its first and
% Second Derivative.
%
%          y = sign(Gamma(x))
%          Gamma   = \sum_{i=1}^{N_sv}\alpha_iy_ik(x,x_i) + b 
%          DGamma  = \sum_{i=1}^{N_sv}-1/2\sigma^2\alpha_iy_ik(x,x_i)(x-x_i)
%          DDGamma = ...
% 
%  The evaluated SVM model can be learned with any toolbox: 
%  libSVM, SMVlight, EnsembleSVM, etc as long as the user creates a 
%  simplified struct model  with the following fields:
%
%  model.nClass: # of Classes (2 for binary)
%  model.nSV    : Total # of Support Vectors
%  model.b      : Offset for classification function
%  model.sigma  : Gaussian RBF kernel Width
%  model.yalphas: Values for the Lagrangian multipliers*Label per SVs [1xnSV]
%  model.SVs    : Set of Support Vectors                              [DxnSV]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load 2D example Dataset and model learned through libSVM
clc; clear all; close all;
load('./mat/2d-example2.mat')

%% Create Simplified Struct Model for SVMGrad from libSVM Model
svmgrad = [];
svmgrad.nClass  = model.nr_class;
svmgrad.nSV     = model.totalSV;
svmgrad.b       = -model.rho;
svmgrad.sigma   = options.sigma;
svmgrad.yalphas = model.sv_coef'; %\alpha_*y_i
svmgrad.SVs     = full(model.SVs)';

%% Visualize Decision Function and gradients (Only for 2d dataset)
plot_svmgrad_boundary(X, labels, svmgrad,  'draw');

%% Sample classifier and gradient evaluation for on query point
query_point = X(1,:)';
tic;
value       = calculateClassifier( svmgrad,  query_point);
gradient    = calculateClassifierDerivative( svmgrad, query_point);
toc;

%% Write SVMGrad Struct to .txt file for C++ Usage
filename = './mat/2d-example.txt';
writeSVMGrad(svmgrad, filename);
