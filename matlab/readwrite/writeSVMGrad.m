function  [] = writeSVMGrad(svmgrad, filename)
%WRITESVMGRAD writes an svmgrad object to a text file
% o svmgrad : svmgrad object
% o filename: filename for text file
%
% The text file will follow the same order of variables
%  model.D       : Datapoint Dimension
%  model.nSV     : Total # of Support Vectors
%  model.b       : Offset for classification function
%  model.sigma   : Gaussian RBF kernel Width
%  model.yalphas : Values for the Lagrangian multipliers*class  [1xnSV]
%  model.SVs     : Set of Support Vectors                       [DxnSV]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(filename,'w');


fclose(fileID);

end

