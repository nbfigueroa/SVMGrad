function [p, tempx, tempy, tempz, classifier_val] = display3DSVM( learned_svm, varargin )
% This function plots the desired components of an SVM object with
% 3-dimensional state space.
%
%   Inputs ----------------------------------------------------------------
%   o learned_svm    :  The SVM object (struct)
%   o varargin   
%       o range        : [xmin xmax ymin ymax zmin zmax] 
%                        'axis' reads the limits of the current axes (gca)
%       o step         : Stepsize for calculating contours
%       o boundary     : 0/1. Select if the classifier boundary is to be drawn
%       o sv           : 0/1. Select if the contours are to be drawn
%       o contour_list : Specific contour value list to be drawn
%       o alpha        : Transparency of the iso-surface patch
%       o targets      : 0/1. Select if the target is to be drawn
%       o color        : Color of the iso-surface patch
%
%   Outputs ---------------------------------------------------------------
%   o p   :  Handle of the 3-D patch object draw as an isosurface
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

ranges=[-1 + learned_svm.target(1) 1 + learned_svm.target(1) -1+ learned_svm.target(2) 1+ learned_svm.target(2) -1+ learned_svm.target(3) 1+ learned_svm.target(3)];
step=0.1;
boundary=false;
targets=true;
svcontour=false;
plotsv=false;
contour_list=[];
alpha=1.0;
colr='g';

for i=1:length(varargin)
    if(strcmp(varargin{i} ,'range'))
        if(strcmp(varargin{i+1}, 'axis'))
            ranges = [get(gca,'XLim'),get(gca,'YLim'),get(gca,'ZLim')];
        else
            ranges=varargin{i+1};
        end
    end
    
    if(strcmp(varargin{i} ,'step'))
        step=varargin{i+1};
    end
    
    if(strcmp(varargin{i} ,'boundary'))
        boundary=varargin{i+1};
    end
    
    if(strcmp(varargin{i} ,'sv_contour'))
        svcontour=varargin{i+1};
    end
    
    if(strcmp(varargin{i} ,'sv'))
        plotsv=varargin{i+1};
    end
    
    if(strcmp(varargin{i} ,'contour_list'))
        contour_list=varargin{i+1};
    end
    
    if(strcmp(varargin{i} ,'alpha'))
        alpha=varargin{i+1};
    end
    if(strcmp(varargin{i} ,'target') || strcmp(varargin{i} ,'targets'))
        targets=varargin{i+1};
    end
    if(strcmp(varargin{i} ,'color'))
        colr=varargin{i+1};
    end
end

    

contour_list = [contour_list,  boundary + 1e-8];


hold on

if(~isempty(contour_list))
    
classifier_val = [];
[ tempx , tempy, tempz ] = meshgrid( ranges(1):step:ranges(2), ranges(3):step:ranges(4),ranges(5):step:ranges(6) );


min_v = inf;
max_v = -inf;
for i = 1 : size ( tempx , 1 )
    for j = 1 : size ( tempx , 2 )
        for k = 1 : size ( tempx , 3 )
            curr_point = [tempx(i,j,k);tempy(i,j,k);tempz(i,j,k)];
            %         curr_point = curr_point/norm(curr_point);
            classifier_val(i,j,k) = calculateClassifier(learned_svm, curr_point);
            if(min_v > classifier_val(i,j,k))
                min_v = classifier_val(i,j,k);
            end
            if(max_v < classifier_val(i,j,k))
                max_v =  classifier_val(i,j,k);
            end
        end
    end
end

    for i=contour_list
        p = patch(isosurface(tempx, tempy, tempz, classifier_val,i));
        set(p,'FaceColor',colr,'EdgeColor','none','FaceAlpha',alpha);
    end
end

if(plotsv)
    for i=1:length(learned_svm.alpha)
        plot3(learned_svm.Sva(1,i), learned_svm.Sva(2,i),learned_svm.Sva(3,i),'ko','Markersize',10,'Linewidth',2);
    end
    for i=1:length(learned_svm.beta)
        plot3(learned_svm.Svb(1,i), learned_svm.Svb(2,i), learned_svm.Svb(3,i),'k^','Markersize',10,'Linewidth',2);
    end
  
end

if(svcontour)
       p = patch(isosurface(tempx, tempy, tempz, classifier_val,1.000001));
%     isonormals(tempx, tempy, tempz, classifier_val, p);
    set(p,'FaceColor',[0 1 0],'EdgeColor','none','FaceAlpha',alpha);
%     daspect([1,1,1])
%     view(3); axis tight

end
% 

if(targets)
plot3(learned_svm.target(1), learned_svm.target(2), learned_svm.target(3), 'ko','Linewidth',2,'Markersize',12,'MarkerFaceColor','m');
end
    camlight
    lighting gouraud
end

