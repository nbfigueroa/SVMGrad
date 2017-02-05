%% 36D Real Robot Collision Data
clear all; close all;
load('Fender-Collision-Avoidance-Dataset.mat')

%% Load CPSP model learnt from this dataset
load('./models/cpsp_models/robot-collision-540k.mat')
clear cpsp_model
cpsp_model = cpsp_model_540k;

%% Make SVM-CPSP model to SVMGrad
svmgrad_cpsp         = [];
svmgrad_cpsp.D       = cpsp_model.D;
svmgrad_cpsp.nSV     = cpsp_model.Nbvs;
svmgrad_cpsp.b       = cpsp_model.bias;
svmgrad_cpsp.sigma   = sqrt(1/(2*cpsp_model.gamma));
svmgrad_cpsp.yalphas = cpsp_model.alphay'; %\alpha_*y_i
svmgrad_cpsp.SVs     = cpsp_model.BVs';

%% Visualize Decision Function and gradients (Only for 2d dataset)
plot_svmgrad_boundary(X, labels, svmgrad_cpsp,  'draw');

%% Sample classifier and gradient evaluation for on query point
query_point = X_test(randi(length(X_test)),:)';
tic;
class       = calculateClass( svmgrad_cpsp,  query_point)
value       = calculateGamma( svmgrad_cpsp,  query_point)
gradient    = calculateGammaDerivative( svmgrad_cpsp, query_point)
toc;

%% Write SVMGrad Struct to .txt file for C++ Usage
filename = './models/Fender/36D-540k-CPSP-Model-Fender.txt';
writeSVMGrad(svmgrad_cpsp, filename);

%% Write Testing Data for SVMGRad
filename = './models/Fender/36D-540k-CPSP-Data-Fender.txt';
ntest    = 500;
randidx  = randperm(length(X_test));
x_test   = X_test(randidx(1:ntest),:)';
y        = zeros(1, ntest);
value    = zeros(1, ntest);
gradient = zeros(svmgrad_cpsp.D, ntest);
for i=1:ntest
    query_point      = x_test(:,i);
    y(1,i)           = calculateClass( svmgrad_cpsp,  query_point);
    value(1,i)       = calculateGamma( svmgrad_cpsp,  query_point);
    gradient(:,i)    = calculateGammaDerivative( svmgrad_cpsp, query_point);
end

writeSVMGradTestData(x_test, y, value, gradient, filename)
