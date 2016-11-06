# This Makefile is used under Linux


# for Mac
# MATLABDIR ?= /opt/local/matlab

CXX = g++
#CXX = g++-4.1
LOCAL_INCLUDE_DIRS=$(shell pwd)/../../include

CFLAGS = -Wall -O3 -fPIC -I$(MATLABDIR)/extern/include -I. -I$(LOCAL_INCLUDE_DIRS)  -Wno-unused-result
LDFLAGS= $(CXXFLAGS)  -Wl,--rpath -Wl,$(shell pwd)/../../lib

MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC\#$(CXX) CXX\#$(CXX) CFLAGS\#"$(CFLAGS)" CXXFLAGS\#"$(CFLAGS)" -O LDOPTIMFLAGS="$(LDFLAGS)"
# comment the following line if you use MATLAB on 32-bit computer
MEX_OPTION += -largeArrayDims -L../../lib -lASVMLearning 
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)
vpath %.h $(LOCAL_INCLUDE_DIRS)

all:	matlab

matlab:	binary

binary:  svmtrain.$(MEX_EXT) 

svmtrain.$(MEX_EXT):       svmtrain.c svm.h svm_model_matlab.o
	$(MEX) $(MEX_OPTION) svmtrain.c svm_model_matlab.o

svm_model_matlab.o:     svm_model_matlab.c svm.h
	$(CXX) $(CFLAGS) -c svm_model_matlab.c


clean:
	rm -f *~ *.o *.mex* *.obj svm.o
