

disp('Building SVMUtil object files...');

mex -c -largeArrayDims -O  -I. -I../../include ../SVMWidget/mxUtil.cpp ../../src/asvm.cpp ../../src/util.cpp 
 
  mex  -largeArrayDims -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I../../include mx_calculate_classifier.cpp  asvm.o util.o mxUtil.o

  mex  -largeArrayDims -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I../../include mx_calculate_classifier_derivative.cpp  asvm.o util.o mxUtil.o

  mex  -largeArrayDims -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I../../include mx_get_bsb_value_bsb_derivative_kernel.cpp  asvm.o util.o 

  mex  -largeArrayDims -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I../../include mx_get_osb_value_osb_derivative_kernel.cpp  asvm.o util.o 

  mex  -largeArrayDims -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I../../include mx_get_osb_value_osb_derivative_kernel_dsigma.cpp  asvm.o util.o 

 mex  -largeArrayDims -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I../../include mx_get_bsb_value_bsb_derivative_kernel_dsigma.cpp  asvm.o util.o