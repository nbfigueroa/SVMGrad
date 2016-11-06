function [ v_modulated ] = getModulatedVelocity( curr_svm, v_nominal, point )
% This function calculates the modulated velocity from the nominal velocity
% and learned SVM object at the given point
%
%   Inputs ----------------------------------------------------------------
%   o curr_svm     :  The SVM object (struct)
%   o v_nominal    :  D x 1 nominal velocity vector
%   o point        :  D x 1 vector for the current point in state space
%
%   Outputs ---------------------------------------------------------------
%   o v_modulated  :  D x 1 modulated velocity vector
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

%% Stabilizing the nominal velocity
stabilization_constant = 0.3*norm(v_nominal);
bad_dir = curr_svm.target - point;
bad_dir = bad_dir/norm(bad_dir);
bad_component = v_nominal'*bad_dir;
tmp = v_nominal - bad_component*bad_dir;
v_nominal = tmp + max(bad_component, stabilization_constant)*bad_dir;



%% Calculating Modulation
grad_h = calculateClassifierDerivative(curr_svm, point);
grad_h = grad_h/(norm(grad_h));
v_perpendicular = v_nominal'*grad_h;
v_parallel = v_nominal - v_perpendicular*grad_h;

MAX_TOL=0.1*norm(v_nominal);
MIN_TOL =0.01*norm(v_nominal);
nrm=norm(v_nominal);
curr_val = calculateClassifier(curr_svm, point);
if(curr_val>0)
    boundary_tolerance = MIN_TOL;
    
else
    
    if(curr_val >-1)
        boundary_tolerance = MAX_TOL;
        %                     boundary_tolerance = (MAX_TOL-MIN_TOL)/2*cos((curr_val+1)*pi/2) + (MAX_TOL+MIN_TOL)/2;
        
    else
        boundary_tolerance = MAX_TOL;
    end
end

v_modulated = v_parallel + max(v_perpendicular, boundary_tolerance)*grad_h;
%             v_modulated = v_parallel + max(v_perpendicular, 0.75)*grad_h;



end

