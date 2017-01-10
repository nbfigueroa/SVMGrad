function  [] = writeSVMGrad(svmgrad, filename)
%WRITESVMGRAD writes an svmgrad object to a text file
% o svmgrad : svmgrad object
% o filename: filename for text file
%
% The text file will follow the same order of variables
%  model.nClass: # of Classes (2 for binary)
%  model.nSV   : Total # of Support Vectors
%  model.b     : Offset for classification function
%  model.sigma : Gaussian RBF kernel Width
%  model.alphas: Values for the Lagrangian multipliers per SVs      [nSV]
%  model.y     : Labels corresponding to SVs                        [nSV]
%  model.SVs   : Set of Support Vectors                             [DxnSV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(filename,'w');


fclose(fileID);

end

