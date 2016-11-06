/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * 
 *  mx_calculate_classifier.cpp
 *
 *  Created on : Jan 22, 2013
 *  Author     : Ashwini Shukla
 *  Email      : ashwini.shukla@epfl.ch
 *  Website    : lasa.epfl.ch
 *
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */

#include "../../include/asvm.h"
#include "../SVMWidget/mxUtil.h"
#ifdef MX_API_VER
#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif
#endif

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	if(nrhs < 2 || !mxIsStruct(prhs[0]) || mxGetN(prhs[1]) == 0)
	{
		printf("Usage: mx_calculate_classifier(svm_struct, data_set)\n");
		return;
	}
	asvm* svm = mxUtil_SVM_mat2c(prhs[0]);
	if(!svm)
		return;

	int dim = mxGetM(prhs[1]);
	int npoints = mxGetN(prhs[1]);
	double* datasetarray = (double*)mxGetPr(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix(npoints, 1, mxREAL);
	double* valptr = mxGetPr(plhs[0]);

	for(int i=0;i<npoints;i++)
		valptr[i] = svm->getclassifiervalue(datasetarray+dim*i);
}

