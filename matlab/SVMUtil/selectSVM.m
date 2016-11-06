function [ curr_val, curr_id, decisive ] = selectSVM( svm_list , point)
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

num_svms = length(svm_list);
tmp =zeros(num_svms,1);
for k=1:num_svms
    tmp(k) = calculateClassifier(svm_list{k}, point);
end

decisive=1;
[curr_val, curr_id] = max(tmp);
% if(curr_val > 0.5)        % in the demonstration region
%     decisive = 1;
% else   % outside demonstrations, but a valid classification region
%         tmp1 = inf;
%         for k=1:num_svms
%             if(tmp1 > norm(svm_list{k}.cg - point))
%                 tmp1 = norm(svm_list{k}.cg - point);
%                 curr_id = k;
%             end
%         end
%         curr_val = tmp(curr_id);
%         decisive = -1;
% end
    
end




