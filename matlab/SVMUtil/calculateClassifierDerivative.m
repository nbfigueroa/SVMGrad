function [ gradient ] = calculateClassifierDerivative( currsvm, point )
% This function calculates the gradient of the classifier function at the
% given point
%
%   Inputs ----------------------------------------------------------------
%   o currsvm  :  The SVM object (struct)
%   o point    :  Vector of length D (dimension of state space)
%
%   Outputs ---------------------------------------------------------------
%   o gradient :  D x 1 gradient vector.
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

gradient =zeros(length(point),1);

d= currsvm.lambda;
type = currsvm.type;
dim = length(point);
e = eye(dim);
target = currsvm.target;

for i=1:length(currsvm.alpha)
    gradient  = gradient + currsvm.alpha(i)*currsvm.y(i)*getKernelFirstDerivative(point, currsvm.Sva(1:dim,i),d, type,1);
end

for i=1:length(currsvm.beta)
    gradient = gradient + currsvm.beta(i)*...
        getKernelSecondDerivative(point, currsvm.Svb(1:dim,i), d, type)*currsvm.Svb(dim+1:2*dim,i);
end

for i=1:length(currsvm.gamma)
    gradient = gradient - currsvm.gamma(i)*getKernelSecondDerivative(point, target, d, type)*e(:,i);
end

end

