function [ value ] = calculateClassifier( currsvm, point )
% This function calculates the value of the classifier function at the
% given point
%
%   Inputs ----------------------------------------------------------------
%   o currsvm :  The SVM object (struct)
%   o point   :  Vector of length D (dimension of state space)
%
%   Outputs ---------------------------------------------------------------
%   o value   :  Scalar value of the function evaluated at this point
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

d= currsvm.lambda;
type = currsvm.type;
dim = length(point);
e = eye(dim);
target = currsvm.target;


value = currsvm.b0;
for i=1:length(currsvm.alpha)
    value  = value + currsvm.alpha(i)*currsvm.y(i)*getKernel(point, currsvm.Sva(1:dim,i),d, type);
end


for i=1:length(currsvm.beta)
    value = value + currsvm.beta(i)*currsvm.Svb(dim+1:2*dim,i)'*...
        getKernelFirstDerivative(point, currsvm.Svb(1:dim,i), d, type, 2);
end

for i=1:length(currsvm.gamma)
    value = value - currsvm.gamma(i)*e(:,i)'*getKernelFirstDerivative(point, target, d, type, 2);
end

end

