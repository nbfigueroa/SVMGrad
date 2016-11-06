function [ new_data ] = resampleAndSmoothData( traj_data, smooth_window, num_data )
% This function smoothes the data and resamples using spline fitting
%
%   Inputs ----------------------------------------------------------------
%   o traj_data     :  N x 1 cell of N trajectories (D x M_i) of D-dimension
%			and M_i points long for i=1...N
%   o smooth_window :  Scalar representing the smoothing window size
%   o num_data      :  Scalar representing the resampled length of trajs.
%
%   Outputs ---------------------------------------------------------------
%   o new_data   :  Resampled and smoothed trajectory N x 1 cell.
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

new_data = cell(length(traj_data),1);
    
   for j=1:length(traj_data)
       for k=1:size(traj_data{j},1)
        traj_data{j}(k,:) = smooth(traj_data{j}(k,:)',smooth_window)';
       end
        
       if(num_data ~= -1)
        new_data{j} = spline(linspace(0,1,length(traj_data{j})), traj_data{j}(:,:), linspace(0,1,num_data));
       else
           new_data = traj_data;
       end

   end



end

