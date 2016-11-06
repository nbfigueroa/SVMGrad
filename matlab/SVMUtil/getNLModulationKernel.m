function [ K, G, G_s, H, H_s, H_ss ] = getNLModulationKernel( input_data, target, labels, d, type )
% This function calculates the kernel matrix for A-SVM
%
%   Inputs ----------------------------------------------------------------
%   o input_data :  D x N matrix of N data points of dimension D
%   o target     :  D x 1 vector representing the target location
%   o labels     :  N x 1 label vector of +1/-1
%   o d          :  Scalar kernel parameter
%   o type       :  String kernel type - 'rbf', 'poly' etc.
%
%   Outputs ---------------------------------------------------------------
%   o K, G, H   :  Block matrices which form the full kernel matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Copyright (c) 2012 Ashwini SHUKLA, LASA, EPFL,          %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The program is free for non-commercial academic use. Please contact the
% author if you are interested in using the software for commercial purposes.
% The software must not be modified or distributed without prior permission
% of the authors. Please acknowledge the authors in any academic publications
% that have made use of this code or part of it. Please use this BibTex
% reference:
%
%
% To get latest upadate of the software please visit
%                          http://asvm.epfl.ch
%
% Please send your feedbacks or questions to:
%                           ashwini.shukla_at_epfl_dot_ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pl = find(labels==1);

M = size(input_data,2);
N = size(input_data,1)/2;
P = length(pl);

K = zeros(M,M);
G = zeros(M,P);
G_s = zeros(M,N);
H = zeros(P,P);
H_s = zeros(P,N);
H_ss = zeros(N,N);

e = eye(N);

for i=1:M
    yi = labels(i);
    xi = input_data(1:N,i);
    for j=1:M
        yj = labels(j);
        xj = input_data(1:N,j);
        K(i,j) = yi*yj*getKernel(xi, xj, d, type);
    end
end

for i=1:M
    yi = labels(i);
    xi = input_data(1:N,i);
    for j=1:P
        xj = input_data(1:N,pl(j));
        xjdot = input_data(N+1:end,pl(j));
        G(i,j) = yi*getKernelFirstDerivative(xi,xj, d, type, 2)'*xjdot;
    end
end


for i=1:M
    xi = input_data(1:N,i);
    yi = labels(i);
    for j=1:N
        G_s(i,j) = yi*getKernelFirstDerivative(xi, target, d, type, 2)'*e(:,j);
    end
end

for i=1:P
    xi = input_data(1:N,pl(i));
    xidot = input_data(N+1:end,pl(i));
    for j=1:P
        xj = input_data(1:N,pl(j));
        xjdot = input_data(N+1:end,pl(j));
        H(i,j) = xidot'*getKernelSecondDerivative(xi, xj, d, type)*xjdot;
    end
end

for i=1:P
    xi = input_data(1:N, pl(i));
    xidot = input_data(N+1:end,pl(i));
    for j=1:N
        H_s(i,j) = xidot'*getKernelSecondDerivative(xi, target, d, type)*e(:,j);
    end
end

for i=1:N
    for j=1:N
        H_ss(i,j) = e(:,i)'*getKernelSecondDerivative(target, target, d, type)*e(:,j);
    end
end

