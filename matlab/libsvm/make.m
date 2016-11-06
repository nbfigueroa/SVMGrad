% To be used if cmake from the root directory fails OR on windows.

disp('Building svmtrain...');
mex   -largeArrayDims -I../../include -O CFLAGS='w' svmtrain.cpp svm_model_matlab.cpp ../../src/svm.cpp
