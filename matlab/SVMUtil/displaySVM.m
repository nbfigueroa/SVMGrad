function [handles, tempx, tempy, classifier_val] =  displaySVM( learned_svm, varargin )
% This function plots the desired components of an SVM object with
% 2-dimensional state space.
%
%   Inputs ----------------------------------------------------------------
%   o learned_svm    :  The SVM object (struct)
%   o varargin
%       o range        : [xmin xmax ymin ymax zmin zmax]
%                        'axis' reads the limits of the current axes (gca)
%       o step         : Stepsize for calculating contours
%       o boundary     : Contour value to be treated as boundary
%       o sv           : 0/1. Select if the contours are to be drawn
%       o num_contours : Number of contours for the positive class
%       o colormap     : Colormap type
%       o targets      : 0/1. Select if the target is to be drawn
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
ranges=[-1 + learned_svm.target(1) 1 + learned_svm.target(1) -1+ learned_svm.target(2) 1+ learned_svm.target(2) ];
step=0.1;
boundary=0;
targets=false;
svcontour=false;
plotsv=false;
num_contours = 10;
clrmp = 'autumn';

for i=1:2:length(varargin)
    if(strcmp(varargin{i} ,'range'))
        if(strcmp(varargin{i+1}, 'axis'))
            ranges = [get(gca,'XLim'),get(gca,'YLim'),get(gca,'ZLim')];
        else
            ranges=varargin{i+1};
        end

        
    else if(strcmp(varargin{i} ,'step'))
            step=varargin{i+1};
            
            
        else if(strcmp(varargin{i} ,'boundary'))
                boundary=varargin{i+1};
                
                
            else if(strcmp(varargin{i} ,'sv_contour'))
                    svcontour=varargin{i+1};
                    
                    
                else if(strcmp(varargin{i} ,'sv'))
                        plotsv=varargin{i+1};
                        
                        
                    else if(strcmp(varargin{i} ,'target') || strcmp(varargin{i} ,'targets'))
                            targets=varargin{i+1};
                            
                        else if(strcmp(varargin{i} ,'colormap') )
                                clrmp=varargin{i+1};
                            else if(strcmp(varargin{i} ,'num_contours'))
                                    num_contours=varargin{i+1};
                                else
                                    disp(['WARNING: Unrecognized property name "' varargin{i} '" !!']);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
hold on

xr = ranges(1) : step : ranges(2);
yr = ranges(3) : step : ranges(4);
if(xr(end) ~= ranges(2))
    xr = [xr ranges(2)];
end
if(yr(end) ~= ranges(4))
    yr = [yr ranges(4)];
end
[ tempx , tempy ] = meshgrid( xr,yr );

hold on
% 
min_v = inf;
max_v = -inf;
classifier_val = zeros(size(tempx));
for i = 1 : size ( tempx , 1 )
    for j = 1 : size ( tempx , 2 )
        curr_point = [tempx(i,j);tempy(i,j)];
        %         curr_point = curr_point/norm(curr_point);
        classifier_val(i,j) = calculateClassifier(learned_svm, curr_point);
        if(min_v > classifier_val(i,j))
            min_v = classifier_val(i,j);
        end
        if(max_v < classifier_val(i,j))
            max_v =  classifier_val(i,j);
        end
    end
end


handle_count = 1;
handles = [];


hl=pcolor(tempx, tempy, classifier_val);
set(get(get(hl,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
shading interp;
colormap(clrmp);
ctr_range = linspace(0,max_v,num_contours);
[C,h]=contour(tempx, tempy, classifier_val,ctr_range(2:end), 'Color','k' );
% [tmp, handles(handle_count)] = contour(tempx, tempy, classifier_val, 1e-10 - boundary,'Linewidth',3,'Color','k','Linestyle','--');

[tmp, handles(handle_count)] = contour(tempx, tempy, classifier_val, boundary, 'Linewidth',3,'Color','k','Linestyle','--');
handle_count = handle_count+1;


if(plotsv)
    hl=[];
    for i=1:length(learned_svm.alpha)
        hl = [hl plot(learned_svm.Sva(1,i), learned_svm.Sva(2,i),'ko','Markersize',12,'Linewidth',2)];
    end
    asv_leg = hggroup;
    set(hl, 'Parent',asv_leg);
    set(get(get(asv_leg,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    hl=[];
    for i=1:length(learned_svm.beta)
        hl = [hl plot(learned_svm.Svb(1,i), learned_svm.Svb(2,i),'k^','Markersize',12,'Linewidth',2)];
    end
    
    bsv_leg = hggroup;
    set(hl, 'Parent',bsv_leg);
    set(get(get(bsv_leg,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    
    handle_count = handle_count+2;
end

if(svcontour)
    min_v = inf;
    for i=1:size(learned_svm.Sva,2)
        if(learned_svm.y(i) ==1)
            tmp = calculateClassifier(learned_svm, learned_svm.Sva(1:2,i));
            if(tmp < min_v)
                min_v = tmp;
            end
        end
        
    end
    [tmp,handles(handle_count)] = contour(tempx, tempy, classifier_val,min_v,'Linewidth',2,'Color','w','LineStyle','-.');
end

if(targets)
    plot(learned_svm.target(1), learned_svm.target(2), 'ko','Markersize',12, 'Linewidth',2, 'MarkerFacecolor','m');
end

