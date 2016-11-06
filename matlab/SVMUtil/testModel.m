function [ traj_list, spurious_att_list, h_val_list ] = testModel( learned_svm, ranges, step, surface_to_mesh, skip_mesh )
% This function tests for spurious attractors in a model by meshing a
% particular iso-surface of the classifier function and simulating trajectories
% starting from the obtained mesh-points.
%
%   Inputs ----------------------------------------------------------------
%   o learned_svm        :  D-dimensional A-SVM model to test (struct)
%   o ranges             :  2D x 1 representing the bounding box for meshing.
%   o step               :  Scalar representing the step size for meshing
%   o surface_to_mesh    :  Scalar representing the iso-surface to mesh.
%   o skip_mesh          :  Scalar representing how to resample mesh-points.
%
%   Outputs ---------------------------------------------------------------
%   o traj_list          :  D x M_i matrix containing all obtained trajectories
%   o spurious_att_list  :  D x N matrix containing N spurious attractors.
%   o h_val_list         :  N x 1 vector containing classifier value at spurious
%			     attractors.
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

target = learned_svm.target;

dim = size(target,1);

spurious_att_list = [];
h_val_list=[];

max_len = 40000;
op = [step,max_len];

switch dim
    case 2

        disp('Meshing surface...');
        [tempx, tempy] = meshgrid(ranges(1):step:ranges(2),ranges(3):step:ranges(4));
        classifier_val = zeros(size(tempx));
        for i = 1 : size ( tempx , 1 )
            for j = 1 : size ( tempx , 2 )
                    curr_point = [tempx(i,j);tempy(i,j)];
                    classifier_val(i,j) = calculateClassifier(learned_svm, curr_point);
            end
        end
        [C,h] = contour(tempx, tempy, classifier_val,[surface_to_mesh surface_to_mesh]);
  
%         set(h,'Color','k');
        vv(:,1) = C(1,2:end)';
        vv(:,2) = C(2,2:end)';
    
        disp(['Found ' num2str(size(vv,1)) ' vertices']);
        disp('Calculating trajectories...');
        u=[];v=[];w=[];
        for i = 1 : size ( tempx , 1 )
            for j = 1 : size ( tempx , 2 )
                    point = [tempx(i,j);tempy(i,j)];
                    
%                     v_nom = learned_svm.target - point;
%                     v_nom = v_nom/norm(v_nom);
%                     direc = getModulatedVelocity(learned_svm, v_nom, point);
                    
                    direc = calculateClassifierDerivative(learned_svm, point);
                    if(norm(direc) >1)
                        direc = direc/norm(direc);
                    end
                    u(i,j) = direc(1);
                    v(i,j) = direc(2);
            end
        end

        
        traj_list = stream2(tempx, tempy, u,v, vv(:,1), vv(:,2),op);
        
        disp('Finding spurious attractors...');
        for i=1:length(traj_list)
            if(~isempty(traj_list{i}) && size(traj_list{i},1) > 4)
                if(norm(traj_list{i}(end,:)' - target) > step && ...
                        abs(traj_list{i}(end,1) - ranges(1)) > step && ...
                        abs(traj_list{i}(end,1) - ranges(2)) > step && ...
                        abs(traj_list{i}(end,2) - ranges(3)) > step && ...
                        abs(traj_list{i}(end,2) - ranges(4)) > step)
                    flag = false;
                    for j=1:size(spurious_att_list,2)
                        if(norm(spurious_att_list(:,j)-traj_list{i}(end,:)') < step)
                            flag = true;
                            break;
                        end
                    end
                    if(~flag)
                        spurious_att_list = [ spurious_att_list, traj_list{i}(end,:)'];
                        h_val_list = [h_val_list;calculateClassifier(learned_svm, point)];
                    end
                end

            end
        end
        disp([num2str(size(spurious_att_list,2)) ' spurious attractors found']);
        
        hold on
        grid on
        
        op=[];
        op.ranges = ranges;
        op.step = step;
        op.plot_boundary = 1;
        %displaySVM(learned_svm,op);
%         for i=1:length(traj_list)
%             if(~isempty(traj_list{i}))
%                 plot(traj_list{i}(:,1), traj_list{i}(:,2),'b-');
%             end
%         end
streamline(traj_list);
        
        if(~isempty(spurious_att_list))
            plot(spurious_att_list(1,:)',spurious_att_list(2,:)','kx','Linewidth',3,'Markersize',18);
        end

        grid on
        axis(ranges);
    case 3
        
        disp(['Meshing surface h = ' num2str(surface_to_mesh)]);
        [tempx, tempy, tempz] = meshgrid(ranges(1):step:ranges(2),ranges(3):step:ranges(4),ranges(5):step:ranges(6));
        classifier_val = zeros(size(tempx));
        
        for i = 1 : size ( tempx , 1 )
            for j = 1 : size ( tempx , 2 )
                for k = 1 : size ( tempx , 3 )
                    curr_point = [tempx(i,j,k);tempy(i,j,k);tempz(i,j,k)];
                    classifier_val(i,j,k) = calculateClassifier(learned_svm, curr_point);
                end
            end
        end
        
%         [tempx, tempy, tempz, classifier_val] = reducevolume(tempx, tempy, tempz, classifier_val, skip_mesh);
        p = patch(isosurface(tempx, tempy, tempz, classifier_val,surface_to_mesh));
        p=reducepatch(p,100);
        vv = p.vertices;
  
        disp(['Found ' num2str(size(vv,1)) ' vertices']);
        disp('Calculating trajectories...');
        u=[];v=[];w=[];
        tt=learned_svm.target;
        
        for i = 1 : size ( tempx , 1 )
            for j = 1 : size ( tempx , 2 )
                for k = 1 : size ( tempx , 3 )
                    curr_point = [tempx(i,j,k);tempy(i,j,k);tempz(i,j,k)];

                    direc = calculateClassifierDerivative(learned_svm, curr_point);
                    if(norm(direc) >0.5)
                        direc = 0.5*direc/norm(direc);
                    end
                    
                    u(i,j,k) = direc(1);
                    v(i,j,k) = direc(2);
                    w(i,j,k) = direc(3);
                end
            end
        end
        
        traj_list = stream3(tempx, tempy, tempz, u,v,w,vv(:,1), vv(:,2), vv(:,3),op);

% dt=step/5;
% traj_list={};
% for i=1:size(vv,1)
%     i
%    curr_point=vv(i,:)';
%    count=1;
%    while(count<100)
%        direc = calculateClassifierDerivative(learned_svm, curr_point(:,end));
%        if(norm(direc) >0.5)
%            direc = 0.5*direc/norm(direc);
%        end
%        curr_point(:,end+1) = curr_point(:,end) + direc*dt;
%        count=count+1;
%        if(norm(curr_point(:,end)) < step)
%           break; 
%        end
%    end
%    traj_list{i} = curr_point';
% end
        
        disp('Finding spurious attractors...');
        for i=1:length(traj_list)
            if(size(traj_list{i},1) >2)
                if(norm(traj_list{i}(end,:)' - target) > step)
                    flag = false;
                    for j=1:size(spurious_att_list,2)
                        if(norm(spurious_att_list(:,j)-traj_list{i}(end,:)') < step)
                            flag = true;
                            break;
                        end
                    end
                    if(~flag)
                        spurious_att_list = [ spurious_att_list, traj_list{i}(end,:)'];
                        h_val_list = [h_val_list;calculateClassifier(learned_svm, traj_list{i}(end,:)')];
                    end
                end

            end
        end
        
        
        disp([num2str(size(spurious_att_list,2)) ' spurious attractors found']);

        
        figure
        hold on
        grid on
        op=[];
        op.ranges = ranges;
        op.step = step;
        op.alpha=0.4;
        op.contour_list = surface_to_mesh;
        display3DSVM(learned_svm,'boundary',false, 'contour_list',surface_to_mesh, 'alpha',0.4, 'range',ranges, 'step',step);
        for i=1:length(traj_list)
            if(~isempty(traj_list{i}))
                plot3(traj_list{i}(:,1), traj_list{i}(:,2),traj_list{i}(:,3),'b-');
                plot3(traj_list{i}(1,1),traj_list{i}(1,2),traj_list{i}(1,3),'ko','Linewidth',2);
            end
        end
        
        if(~isempty(spurious_att_list))
            plot3(spurious_att_list(1,:)',spurious_att_list(2,:)',spurious_att_list(3,:)','kx','Linewidth',3,'Markersize',18);
        end
        
        axis(ranges);
    case default
        
end

end


